# Loading needed libraries

import warnings

import numpy as np
np.float_ = np.float64
import pandas as pd

with warnings.catch_warnings(): # Filtering harmless warnings
    warnings.filterwarnings(
        "ignore",
        message=".*Dask DataFrame.*"
    )
    warnings.filterwarnings(
        "ignore",
        message=".*Import anndata.io.read_text instead.*"
    )
    import squidpy as sq

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, to_hex

from itertools import combinations, product
import networkx as nx
from collections import defaultdict, deque
import logging

logger = logging.getLogger(__name__)

def _expand_subsets(col,df):
    """
    Expands the subsets of a cell type based on the presence or absence of protein markers.

    Parameters
    -----
    col : str
        Column name in the DataFrame.
    df : pd.DataFrame
        DataFrame containing the cell type information.

    Returns
    ------
    pd.DataFrame
       Expanded DataFrame with subsets of the cell type.
    """
    n = df[col].apply(lambda x: isinstance(x, list)).sum()
    combos = list(product(*([ [0,1] ] * n)))
    list_to_df = []
    aux = 0
    for val in df[col]:
        if isinstance(val,list):
            list_to_df.append([c[aux] for c in combos])
            aux+=1
        else: list_to_df.append([val]*2**n)
    # return pd.DataFrame(list_to_df, columns=[f"{col}_{''.join(map(str, c))}" for c in combos])
    return pd.DataFrame(list_to_df, columns=[f"{col}_{c}" for c in range(len(combos))])

def _move_bigger_groups(dic, zero_degree_items, target, final_subsets):
    """
    Reorganize less-specific phenotype groups during subset resolution.
    """
    find = []
    for k, item_list in dic.items():
        for item in zero_degree_items:
            if item in item_list:
                item_list.remove(item)
                find.append(item)
    dic[target].extend(find)
    
    # Number of None values will determine the order of the cell types. The more None values, the bigger the group.
    dicc_target = {k:[] for k in reversed(range(target+1))}
    dicc_target[target] = dic[target]
    for i in reversed(range(target)):
        # if dic[1]!=0:
        for item in dic[i]:
            check = 0
            for k, item_list in final_subsets.items():
                if item in item_list:
                    for keys, vals in dicc_target.items():
                        if k in vals:
                            check+=1
                            dicc_target[keys-1].append(item)
                            break
            if check == 0: dicc_target[target].append(item)
    for keys, vals in dicc_target.items():
        dicc_target[keys] = list(dict.fromkeys(vals)) # Remove duplicates while preserving order
        for val in vals:
            if val in final_subsets.keys():
                dicc_target[keys-1].extend(final_subsets[val])
    dicc_target = {k:v for k,v in dicc_target.items() if len(v)!=0}
    if len(dicc_target)<target+1:
        for i in range(target+1):
            if i not in dicc_target:
                for j in range(i+1, target+1):
                    if j in dicc_target:
                        dicc_target[i] = dicc_target[j]
                        dicc_target.pop(j)
                        break
    return dicc_target, max([k for k in dicc_target.keys()])

def _mix_hex_colors(colors,dim_factor=0.25,background="#ffffff"):
    """
    Combine multiple HEX colors and optionally dim/desaturate the result
    by blending it toward a background color.

    Parameters
    ----------
    colors : list[str]
        HEX colors to combine.
    dim_factor : float
        Amount of blending toward the background.
        0.0 keeps the original mixed color.
        1.0 becomes the background color.
    background : str
        Background color used for dimming.
        Use '#ffffff' for white plots or '#000000' for dark plots.

    Returns
    -------
    str
        Combined HEX color.
    """
    rgb_colors = [to_rgb(color) for color in colors]
    mixed_rgb = tuple(sum(color[channel] for color in rgb_colors) / len(rgb_colors) for channel in range(3))
    background_rgb = to_rgb(background)
    dimmed_rgb = tuple((1 - dim_factor) * mixed + dim_factor * background_value for mixed, background_value in zip(mixed_rgb, background_rgb))
    return to_hex(dimmed_rgb)

def _create_intersection_color(intersection_label, list_intersected_cell_subtypes, custom_colors, dim_factor=0.20, background="#ffffff"):
    """
    Generate a color for an intersection label from its component colors.

    Parameters
    ----------
    intersection_label : str
        Cell type that comes from the intersection of two cell types.
    list_intersected_cell_subtypes: List[str]
        List with the cell types who are intersected.
    custom_colors : Dict[str:str]
       Dictionary containing HEX colors per cell type. 
    dim_factor : float
        Amount of blending toward the background.
        0.0 keeps the original mixed color.
        1.0 becomes the background color.
    background : str
        Background color used for dimming.
        Use '#ffffff' for white plots or '#000000' for dark plots.

    Returns
    -------
    str
        Combined HEX color.
    """
    if intersection_label not in custom_colors.keys():
        component_colors = [custom_colors[subtype] for subtype in list_intersected_cell_subtypes]
        return _mix_hex_colors(component_colors, dim_factor=dim_factor, background=background)
    else: return custom_colors[intersection_label]

def plot_digraph(ct_df, protein_markers, output_dir, output_filename_identifier, order_dicc, custom_colors):
    """
    Resolve phenotype subset relationships and plot their directed graph.

    Parameters
    ----------
    ct_df : pandas.DataFrame
        Marker-rule table whose columns represent phenotypes.
    protein_markers : list of str
        Marker names represented in the rule table.
    output_dir : pathlib.Path
        Base SpPrAn output directory.
    output_filename_identifier : str
        Label used to identify the generated hierarchy output.
    order_dicc : dict
        Phenotype ordering grouped by rule specificity.
    custom_colors : dict of str to str
        Phenotype color mapping.

    Returns
    -------
    final_subsets : dict
        Resolved subset relationships.
    subtype_definitions : dict
        Phenotype definitions generated during subset resolution.
    larger_groups : list or dict
        Less-specific parent groups retained by the hierarchy.
    custom_colors : dict
        Updated phenotype color mapping.

    Notes
    -----
    The directed graph represents logical set containment among
    marker-defined phenotypes rather than a biological developmental lineage.
    """
    combinations_list = [list(c) for c in combinations(ct_df.columns, 2)]
    dicc_subsets = {col:_expand_subsets(col, ct_df) for col in ct_df.columns}

    dicc_temp = {}
    max_N = max([k for k in order_dicc.keys()])
    for total_N, cell_types in order_dicc.items():
        dicc_temp[total_N] = []
        for cell_type in cell_types:
            if cell_type in ct_df.columns:
                dicc_temp[total_N].append(cell_type)
        if total_N==max_N and len(dicc_temp[total_N])==0:
            dicc_temp.pop(total_N)
            max_N = max([k for k in order_dicc.keys() if k<total_N])

    cell_type_chart = nx.DiGraph()
    for cell_type in ct_df.columns: cell_type_chart.add_node(cell_type) 
    
    final_subsets = {cell_type:[] for cell_type in ct_df.columns}
    final_subsets_mirror = {cell_type:[] for cell_type in ct_df.columns}
    for elem in combinations_list:
        temp_0 = list(dicc_subsets[elem[0]].columns)
        temp_1 = [col for col in dicc_subsets[elem[0]].columns for col_c in dicc_subsets[elem[1]].columns if dicc_subsets[elem[0]][col].equals(dicc_subsets[elem[1]][col_c])]
        match len(temp_1):
            case x if x==len(temp_0): 
                logger.warning(f'{elem[0]} is subset of {elem[1]}')
                cell_type_chart.add_edges_from([(elem[0],elem[1])])
                final_subsets[elem[1]].append(elem[0])
                if len(temp_1) != 1:
                    for pos_subset in range([key for key in dicc_temp.keys() if elem[0] in dicc_temp[key]][0]):
                        for pot_subset in dicc_temp[pos_subset]:
                            if nx.has_path(cell_type_chart,source=pot_subset,target=elem[1]) and nx.shortest_path_length(cell_type_chart, source=pot_subset, target=elem[1])==1:
                                cell_type_chart.remove_edge(pot_subset, elem[1])
                                final_subsets[elem[1]].remove(pot_subset)
            case x if x>0: 
                whos_who = []
                for pos_subset in range([key for key in dicc_temp.keys() if elem[0] in dicc_temp[key]][0]):
                    for pot_subset in dicc_temp[pos_subset]:
                        if nx.has_path(cell_type_chart,source=pot_subset,target=elem[0]) and nx.has_path(cell_type_chart,source=pot_subset,target=elem[1]):
                            whos_who.append(pot_subset)
                list_temp_temp = list(set(col for col in temp_1 for col_s in whos_who for col_c in dicc_subsets[col_s].columns if dicc_subsets[elem[0]][col].equals(dicc_subsets[col_s][col_c])))
                temp_new = [x for x in temp_1 if x not in list_temp_temp]
                limit_dim_color = {k:0.2+0.05*k for k in range(len(temp_new))} if len(temp_new) <= 14 else {k:0.2+(0.7*k/len(temp_new)) for k in range(len(temp_new))}
                0.2+len(temp_new)
                match len(temp_new):
                    case x if x>1:
                        logger.warning(f"There are sub cell types that are not defined by user and belong to both {elem[0]} and {elem[1]}.")
                        
                        for num_to_add,sct_def in enumerate(temp_new):
                            if sct_def in cell_type_chart:
                                if not nx.has_path(cell_type_chart,source=sct_def,target=elem[0]): 
                                    cell_type_chart.add_edges_from([(sct_def,elem[0])])
                                    final_subsets[elem[0]].append(sct_def)
                                    name_to_add = f"Rare {elem[0]} & {elem[1]}_{num_to_add+1}"
                                    final_subsets_mirror[elem[0]].append((sct_def,name_to_add))
                                if not nx.has_path(cell_type_chart,source=sct_def,target=elem[1]): 
                                    cell_type_chart.add_edges_from([(sct_def,elem[1])])
                                    final_subsets[elem[1]].append(sct_def)
                                    name_to_add = f"Rare {elem[0]} & {elem[1]}_{num_to_add+1}"
                                    final_subsets_mirror[elem[1]].append((sct_def,name_to_add))
                            else: 
                                cell_type_chart.add_node(sct_def)
                                cell_type_chart.add_edges_from([(sct_def,elem[0])])
                                cell_type_chart.add_edges_from([(sct_def,elem[1])])
                                final_subsets[elem[0]].append(sct_def)
                                final_subsets[elem[1]].append(sct_def)
                                name_to_add = f"Rare {elem[0]} & {elem[1]}_{num_to_add+1}"
                                final_subsets_mirror[elem[0]].append((sct_def,name_to_add))
                                final_subsets_mirror[elem[1]].append((sct_def,name_to_add))
                            custom_colors[f"Rare {elem[0]} & {elem[1]}_{num_to_add+1}"] = _create_intersection_color(f"Rare {elem[0]} & {elem[1]}_{num_to_add+1}",[elem[0],elem[1]],custom_colors,limit_dim_color[temp_new.index(sct_def)])
                    case x if x==1:
                        logger.warning(f"There is a sub cell type that is not defined by user and belongs to both {elem[0]} and {elem[1]}.")
                        if temp_new[0] in cell_type_chart:
                            if not nx.has_path(cell_type_chart,source=temp_new[0],target=elem[0]): 
                                cell_type_chart.add_edges_from([(temp_new[0],elem[0])])
                                final_subsets[elem[0]].append(temp_new[0])
                                name_to_add = f"Rare {elem[0]} & {elem[1]}_1"
                                final_subsets_mirror[elem[0]].append((temp_new[0],name_to_add))
                            if not nx.has_path(cell_type_chart,source=temp_new[0],target=elem[1]): 
                                cell_type_chart.add_edges_from([(temp_new[0],elem[1])])
                                final_subsets[elem[1]].append(temp_new[0])
                                name_to_add = f"Rare {elem[0]} & {elem[1]}_1"
                                final_subsets_mirror[elem[1]].append((temp_new[0],name_to_add))
                        else:
                            cell_type_chart.add_node(temp_new[0])
                            cell_type_chart.add_edges_from([(temp_new[0],elem[0])])
                            cell_type_chart.add_edges_from([(temp_new[0],elem[1])])
                            final_subsets[elem[0]].append(temp_new[0])
                            final_subsets[elem[1]].append(temp_new[0])
                            name_to_add = f"Rare {elem[0]} & {elem[1]}_1"
                            final_subsets_mirror[elem[0]].append((temp_new[0],name_to_add))
                            final_subsets_mirror[elem[1]].append((temp_new[0],name_to_add))
                        custom_colors[f"Rare {elem[0]} & {elem[1]}_1"] = _create_intersection_color(f"Rare {elem[0]} & {elem[1]}_1",[elem[0],elem[1]],custom_colors)
                    case _: pass
            case _: pass
    
    possible_duplicates = {cell_type:list_val for cell_type, list_val in final_subsets.items() if len(list_val)>1}
    for cell_type, list_val in possible_duplicates.items():
        for ct in list_val:
            for ct_temp in [k for k in possible_duplicates.keys() if k!=cell_type]:
                if set([cell_type,ct]).issubset(set(possible_duplicates[ct_temp])):
                    final_subsets[ct_temp].remove(ct)
                    cell_type_chart.remove_edge(ct, ct_temp)
    csv_path = output_dir / "results" / output_filename_identifier
    if not csv_path.exists(): csv_path.mkdir(parents=True, exist_ok=True)
    temu_dicc = {cell_type:[] for cell_type in final_subsets.keys()}
    for cell_type, list_val in final_subsets.items():
        for elem_prime in list_val:
            for elem_to_check in final_subsets_mirror[cell_type]:
                if elem_prime == elem_to_check[0]:
                    temu_dicc[cell_type].append(elem_to_check[1])
    pd.DataFrame({k:[v] for k,v in temu_dicc.items()},index=['Subsets']).T.to_csv(csv_path / "Cell type subsets.csv", index=True)
    
    ct_final_subsets = ct_df.copy()
    final_order = []
    new_columns = {}
    past_keys = []
    cols_ct_final_subsets = list(ct_final_subsets.columns)
    for k,val in final_subsets.items():
        keys_new_columns = list(new_columns.keys())
        columns_to_check = cols_ct_final_subsets+keys_new_columns
        if val:
            for ct in val:
                if ct not in columns_to_check and ct not in past_keys:
                    for elem_to_check in final_subsets_mirror[k]:
                        if ct == elem_to_check[0]:
                            past_keys.append(ct)
                            new_columns[elem_to_check[1]]=dicc_subsets[k][ct].rename(index={idx:pmarker for idx,pmarker in enumerate(protein_markers)})
                            final_order.append(elem_to_check[1])
        final_order.append(k)
    if new_columns: ct_final_subsets = pd.concat([ct_final_subsets,pd.DataFrame(new_columns)],axis=1)
    ct_final_subsets = ct_final_subsets[final_order]
    ct_final_subsets = ct_final_subsets.to_dict()
    df_df_temp = pd.DataFrame(ct_final_subsets, dtype=object)
    none_mask_temp = df_df_temp.isna()
    for col in df_df_temp.columns:
        if none_mask_temp[col].any():
            df_df_temp.loc[df_df_temp[col].isna(), col] = [[0,1] for _ in range(df_df_temp[col].isna().sum())]
    df_df_temp.to_csv(csv_path / "All cell type definitions.csv", index=True)
    
    cell_type_chart = cell_type_chart.subgraph(list(ct_df.columns)).copy() # This line avoids to plot not-defined cell types
    zero_degree_ct = [cell_type for cell_type in ct_df.columns if cell_type_chart.degree[cell_type] == 0]
    if max_N != 0: dicc_temp, max_N = _move_bigger_groups(dicc_temp, zero_degree_ct, max_N, final_subsets)

    return temu_dicc, ct_final_subsets, dicc_temp[max_N], custom_colors

def _plot_cell_hierarchy(
        graph: nx.DiGraph,
        custom_colors,
        save_path,
        *,
        figsize: tuple[float, float] | None = None,
        show_marker_rules: bool = True,
        rules_on_edges: bool = False,
        horizontal_spacing: float = 3.0,
        vertical_spacing: float = 2.2,
        title: str = "Cell phenotype hierarchy",
        dpi: int = 300,
        arrows_orientation: str = "contains",
    ):
    """
    Plot a cell hierarchy whose arrows point from child to parent.

    By default, marker rules are placed inside each node. Set
    rules_on_edges=True to place subtype-defining rules on the edges instead.
    """

    def _forest_layout(
            graph: nx.DiGraph,
            *,
            horizontal_spacing: float = 3.0,
            vertical_spacing: float = 2.2,
            root_gap: float = 1.5,
        ) -> dict[str, tuple[float, float]]:
        """Place every parent above the horizontal center of its descendants."""

        parent_to_child = graph.reverse(copy=False)
        roots = graph.graph.get("root_order", [node for node in graph if graph.out_degree(node) == 0])

        positions: dict[str, tuple[float, float]] = {}
        cursor = 0.0

        def _place_subtree(node: str) -> float:
            nonlocal cursor
            children = list(parent_to_child.successors(node))
            if children:
                child_x = [_place_subtree(child) for child in children]
                x_position = sum(child_x) / len(child_x)
            else:
                x_position = cursor
                cursor += horizontal_spacing
            level = graph.nodes[node]["level"]
            positions[node] = (x_position, level * vertical_spacing)
            return x_position

        for index, root in enumerate(roots):
            _place_subtree(root)
            if index < len(roots) - 1:
                cursor += root_gap * horizontal_spacing

        return positions

    def _wrap_rule_text(rule_text: str, max_items_per_line: int = 2) -> str:
        """Wrap semicolon-separated marker rules across multiple lines."""
        if not rule_text:
            return ""

        items = [item.strip() for item in rule_text.split(";") if item.strip()]
        lines = [
            "\n".join(items[index : index + max_items_per_line])
            for index in range(0, len(items), max_items_per_line)
        ]
        return "\n".join(lines)

    def _short_marker_name(marker: str) -> str:
        """Convert a long positivity-column name into a compact marker name."""
        name = str(marker).strip()
        name = name.replace("Positivity - ", "").replace("* (MV - CYTO)"," (Cytoplasm)").replace("* (MV - NUC)"," (Nucleus)").replace("* (MV Cell)"," (Cell)")
        name = name.replace("(Cytoplasm)","").replace(" (Nucleus)","").replace("(Cell)","")
        return name

    def _format_marker_rules(
            rules,
            *,
            separator: str = "; ",
            omit_null: bool = True,
        ) -> str:
        """Format marker states as compact labels."""

        labels: list[str] = []
        for marker, value in rules.items():
            if value is None and omit_null:
                continue
            marker_name = _short_marker_name(marker)
            if value in (1, True): state = "+"
            elif value in (0, False): state = "−"
            elif value is None: state = "any"
            else: state = f"={value}"
            labels.append(f"{marker_name}{state}")
        return separator.join(labels)

    def _contrasting_text_color(background: str) -> str:
        """Select black or white text according to node-background luminance."""
        try: red, green, blue = to_rgb(background)
        except ValueError: return "black"
        luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
        return "black" if luminance > 0.55 else "white"
    
    positions = _forest_layout(graph, horizontal_spacing=horizontal_spacing, vertical_spacing=vertical_spacing,)
    custom_colors = graph.graph.get("custom_colors", {})
    parent_to_child = graph.reverse(copy=False)
    leaves = [node for node in graph if parent_to_child.out_degree(node) == 0]
    max_depth = max(abs(attributes["level"]) for _, attributes in graph.nodes(data=True))
    if figsize is None: figsize = (max(16.0, 1.8 * len(leaves)), max(7.0, 2.6 * (max_depth + 1)),)
    if arrows_orientation == "contains" : arrowstyle="<|-"
    else: arrowstyle="-|>"
    fig, ax = plt.subplots(figsize=figsize)
    nx.draw_networkx_edges(
        graph,
        positions,
        ax=ax,
        arrows=True,
        arrowstyle=arrowstyle,
        arrowsize=20,
        width=1.4,
        edge_color="#666666",
        node_size=8500,
        min_source_margin=18,
        min_target_margin=18,
    )

    if show_marker_rules and rules_on_edges:
        edge_labels = {(child, parent): _wrap_rule_text(_format_marker_rules(attributes.get("rules", {})), max_items_per_line=1) for child, parent, attributes in graph.edges(data=True)}
        edge_labels = {edge: label for edge, label in edge_labels.items() if label}
        nx.draw_networkx_edge_labels(
            graph,
            positions,
            edge_labels=edge_labels,
            ax=ax,
            rotate=False,
            font_size=7.5,
            label_pos=0.52,
            bbox={
                "boxstyle": "round,pad=0.15",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.9,
            },
        )

    for node, (x_position, y_position) in positions.items():
        attributes = graph.nodes[node]
        node_color = custom_colors.get(node, "#D9E8F5")
        is_root = attributes.get("kind") == "cell_type"

        node_label = str(node)
        if show_marker_rules and not rules_on_edges:
            rule_text = _wrap_rule_text(_format_marker_rules(attributes.get("rules", {})), max_items_per_line=2)
            if rule_text: node_label = f"{node_label}\n{rule_text}"
        ax.text(
            x_position,
            y_position,
            node_label,
            ha="center",
            va="center",
            fontsize=8.5,
            fontweight="bold" if is_root else "normal",
            color=_contrasting_text_color(node_color),
            linespacing=1.25,
            bbox={
                "boxstyle": "round,pad=0.45",
                "facecolor": node_color,
                "edgecolor": "black",
                "linewidth": 2.0 if is_root else 1.0,
            },
            zorder=3,
        )

    levels = sorted({attributes["level"] for _, attributes in graph.nodes(data=True)}, reverse=True)
    ax.set_yticks([level * vertical_spacing for level in levels])
    ax.set_yticklabels(["Level 0: cell type" if level == 0 else f"Level {level}: subtype generation {abs(level)}" for level in levels])
    ax.set_ylabel("Phenotype hierarchy")
    ax.set_title(title, fontsize=15, pad=20)
    ax.grid(axis="y", linestyle=":", alpha=0.35)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.margins(x=0.03, y=0.12)
    ax.set_axisbelow(True)
    for spine in ax.spines.values(): spine.set_visible(False)
    fig.tight_layout()
    if save_path is not None:
        if not save_path.exists(): save_path.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path / "Cell phenotype hierarchy.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def _build_cell_hierarchy_graph(protein_markers, cell_types, cell_sub_types,custom_colors):
    """
    Build a directed cell-phenotype hierarchy from the configuration data.
    """

    graph = nx.DiGraph()
    graph.graph["custom_colors"] = custom_colors
    graph.graph["protein_markers"] = protein_markers
    graph.graph["root_order"] = list(cell_types.keys())

    # The cell-type vectors correspond positionally to protein_markers.
    for cell_type, signature in cell_types.items(): graph.add_node(cell_type, kind="cell_type", rules=signature)

    # Each top-level key is a parent, and every nested key is its child.
    for parent, children in cell_sub_types.items():
        if parent not in graph: graph.add_node(parent, kind="subtype", rules={})
        for child, rules in children.items():
            rules = {} if rules is None else rules
            if child not in graph:
                graph.add_node(child, kind="subtype", rules=dict(rules))
            else:
                # A subtype may already exist because it was previously added
                # as the parent of another subtype generation.
                graph.nodes[child]["kind"] = graph.nodes[child].get("kind", "subtype")
                graph.nodes[child]["rules"] = dict(rules)
            graph.add_edge(child, parent, rules=dict(rules))

    if not nx.is_directed_acyclic_graph(graph):
        cycle = nx.find_cycle(graph)
        raise ValueError(f"The hierarchy contains a cycle: {cycle}")

    roots = list(cell_types.keys())

    # Reverse only for traversal: parent -> children.
    parent_to_child = graph.reverse(copy=False)
    levels: dict[str, int] = {}

    for root in roots:
        distances = nx.single_source_shortest_path_length(parent_to_child, root)
        for node, distance in distances.items():
            proposed_level = -distance
            levels[node] = proposed_level

    nx.set_node_attributes(graph, levels, "level")
    return graph

def plot_general_cell_hierarchy(protein_markers, cell_types, cell_sub_types, custom_colors, save_path, arrows_orientation="contains"):
    """
    Plot the user-defined primary and subtype phenotype hierarchy.

    Parameters
    ----------
    protein_markers : list of str
        Primary protein-marker positivity columns.
    cell_types : dict
        Primary phenotype definitions.
    cell_sub_types : dict
        Parent-child subtype definitions.
    custom_colors : dict of str to str
        Phenotype color mapping.
    save_path : pathlib.Path
        Directory where the hierarchy figure is saved.
    arrows_orientation : {"contains", "isContained"}, default="contains"
        Direction used for hierarchy arrows.

    Returns
    -------
    None
        The hierarchy figure is saved to disk.

    Notes
    -----
    Nodes represent marker-defined populations. Parent-child edges describe
    nested phenotype definitions supplied by the user.
    """
    graph = _build_cell_hierarchy_graph(protein_markers, cell_types, cell_sub_types, custom_colors)
    _plot_cell_hierarchy(graph, custom_colors, save_path, arrows_orientation=arrows_orientation)

def order_sub_cell_types(cell_subtypes):
    """
    Order subtype parents so ancestors are processed before descendants.

    Parameters
    ----------
    cell_subtypes : dict
        User-defined subtype hierarchy.

    Returns
    -------
    list of str
        Parent phenotype names in dependency order. Original insertion order is
        preserved when multiple nodes have the same dependency level.

    Raises
    ------
    ValueError
        If the dependency structure cannot be topologically ordered.
    """
    top_level_keys = list(cell_subtypes.keys())
    top_level_set = set(top_level_keys)
    children_by_parent = defaultdict(list)
    indegree = {key: 0 for key in top_level_keys}

    for parent, child_dict in cell_subtypes.items():
        if not isinstance(child_dict, dict): continue
        for child in child_dict.keys():
            if child in top_level_set:
                children_by_parent[parent].append(child)
                indegree[child] += 1

    queue = deque([key for key in top_level_keys if indegree[key] == 0])
    ordered = []
    while queue:
        parent = queue.popleft()
        ordered.append(parent)
        for child in children_by_parent[parent]:
            indegree[child] -= 1
            if indegree[child] == 0: queue.append(child)

    if len(ordered) != len(top_level_keys):
        remaining = [key for key in top_level_keys if indegree[key] > 0]
        raise ValueError(
            "Cycle detected among cell_subtypes dependencies. "
            f"These keys could not be ordered: {remaining}"
        )
    return ordered

def export_legend_png(ax,output_path,handles, labels, horizontal=False, dpi=300):
    old_legend     = ax.get_legend()
    title          = None
    frameon        = False
    fontsize       = 10
    title_fontsize = False
    marker_scale   = 1.0
    labelspacing   = 0.6
    columnspacing  = 1.2
    handletextpad  = 0.4
    borderpad      = 0.2
    transparent    = True
    if old_legend is not None and old_legend.get_title() is not None:
        existing_title = old_legend.get_title().get_text()
        title = existing_title if existing_title else None
    if horizontal:
        ncol = len(labels)
        fig_width = max(2.5, min(18.0, 0.9 * len(labels)))
        fig_height = 0.75 if title is None else 1.0
    else:
        ncol = 1
        fig_width = max(2.8,(1/72) * max([len(lab) for lab in labels]))
        fig_height = max(1.2, 0.35 * len(labels) + (0.35 if title else 0.0))
    loc = "center"
    fig = plt.figure(figsize=(fig_width, fig_height))
    legend = fig.legend(handles, labels, loc=loc, ncol=ncol, title=title, frameon=frameon, fontsize=fontsize, title_fontsize=title_fontsize, markerscale=marker_scale, labelspacing=labelspacing, columnspacing=columnspacing, handletextpad=handletextpad, borderpad=borderpad,)
    fig.canvas.draw()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", bbox_extra_artists=[legend], pad_inches=0.05, transparent=transparent,)
    plt.close(fig)
