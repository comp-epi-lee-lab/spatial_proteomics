"""Configuration loading and validation utilities for SpPrAn.

This module is responsible for:

1. Safely parsing YAML configuration files.
2. Rejecting duplicate YAML mapping keys.
3. Applying configuration defaults.
4. Validating configuration structure and cross-field rules.
5. Converting workspace paths to ``pathlib.Path`` objects.
6. Verifying that user-defined input and output directories exist.

The public entry point is :func:`load_config`.
"""

import copy
import logging
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import yaml


logger = logging.getLogger(__name__)

REVIEW_YAML_MESSAGE = "Please review the YAML configuration file."

HEX_COLOR_RE = re.compile(r"^#[0-9A-Fa-f]{6}$")

WORKSPACE_FILETYPES = {"tsv", "csv"}

PLOT_OPTIONS = {"none", "all", "highlight", "overview", "custom",}

CUSTOM_PLOT_TYPES = {"none", "all", "highlight", "overview",}

# Add other supported values here if SPPRAN accepts more orientations.
ARROWS_ORIENTATIONS = {"contains", "isContained",}


class ConfigValidationError(ValueError):
    """
    Report one or more invalid SpPrAn configuration settings.

    Parameters
    ----------
    errors : list of str
        Human-readable validation messages collected while checking the YAML
        configuration.

    Attributes
    ----------
    errors : list of str
        Individual configuration errors.
    """

    def __init__(self, errors: list[str]):
        self.errors = errors
        message = ("Invalid YAML configuration:\n" + "\n".join(f"  - {error}" for error in errors) + f"\n\n{REVIEW_YAML_MESSAGE}")
        super().__init__(message)


class UniqueKeyLoader(yaml.SafeLoader):
    """Safe YAML loader that rejects duplicate mapping keys."""


def _construct_unique_mapping(loader: yaml.SafeLoader, node: yaml.MappingNode, deep: bool = False):
    """Construct a YAML mapping while rejecting duplicate keys."""
    mapping = {}

    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)

        try:
            already_present = key in mapping
        except TypeError as exc:
            raise yaml.constructor.ConstructorError("while constructing a mapping", node.start_mark, "found an unhashable mapping key", key_node.start_mark,) from exc

        if already_present:
            raise yaml.constructor.ConstructorError("while constructing a mapping", node.start_mark, f"found duplicate key {key!r}", key_node.start_mark,)

        mapping[key] = loader.construct_object(value_node, deep=deep)

    return mapping


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_unique_mapping,)


def _is_binary_or_null(value: Any) -> bool:
    """Return True only for integer 0, integer 1, or None."""
    # bool must be excluded explicitly because bool is a subclass of int.
    return value is None or (type(value) is int and value in (0, 1))


def _is_binary(value: Any) -> bool:
    """Return True only for integer 0 or integer 1."""
    return type(value) is int and value in (0, 1)


def _valid_marker_name(marker: Any) -> bool:
    """Check the basic expected format of a protein-marker name."""
    return (
        isinstance(marker, str)
        and marker.strip().startswith("Positivity - ")
        and "(MV" in marker
        and marker.strip().endswith(")")
    )


def _check_yaml(config: Any) -> dict:
    """Validate and normalize an already-parsed SPPRAN configuration.

    Validation errors are collected so users can fix several independent
    problems in one pass rather than discovering one error per run.

    Parameters
    ----------
    config : Any
        Object produced by YAML parsing.

    Returns
    -------
    dict
        A validated, normalized copy of the configuration.

    Raises
    ------
    ConfigValidationError
        If one or more settings are invalid.
    """
    if not isinstance(config, dict) or not config: raise ConfigValidationError(["Configuration file is empty or does not contain a YAML dictionary."])

    config = copy.deepcopy(config)
    errors: list[str] = []

    def error(message: str) -> None:
        errors.append(message)

    def warning(message: str) -> None:
        logger.warning("WARNING: %s", message)

    # ========== workspace ==========
    workspace = config.get("workspace")

    if workspace is None:
        error("'workspace' setting is missing.")
        workspace = {}
        config["workspace"] = workspace
    elif not isinstance(workspace, dict):
        error("'workspace' must be a dictionary, " f"got {type(workspace).__name__}.")
        workspace = {}
        config["workspace"] = workspace

    for path_key, label in (("input_dir", "Input"), ("output_dir", "Output")):
        if path_key not in workspace:
            error(f"'workspace.{path_key}' setting is missing.")
            continue

        path_value = workspace[path_key]
        if not isinstance(path_value, str) or not path_value.strip(): error(f"'{label} directory' must contain a non-empty path string.")

    input_dir = workspace.get("input_dir")
    input_dir_checked = None
    output_dir = workspace.get("output_dir")

    if (isinstance(input_dir, str) and input_dir.strip() and isinstance(output_dir, str) and output_dir.strip() and os.path.normpath(input_dir) == os.path.normpath(output_dir)):
        warning("Input and output directories are the same. Using a different output directory is recommended to avoid mixing input and results.")

    if not Path(input_dir).exists(): error(f"Input directory '{input_dir}' does not exist or is misspelled.")
    if not Path(input_dir).is_dir(): error(f"Input path '{input_dir}' exists, but it is not a directory.")
    else: input_dir_checked = Path(input_dir)
    if not Path(output_dir).exists(): error(f"Output directory '{output_dir}' does not exist or is misspelled. Please create the desired output directory and specify it in the configuration file.")
    if not Path(output_dir).is_dir(): error(f"Output path '{output_dir}' exists, but it is not a directory.")

    if "filetype" not in workspace:
        warning("'workspace.filetype' is missing. Defaulting to 'tsv'.")
        workspace["filetype"] = "tsv"
    elif not isinstance(workspace["filetype"], str): error(f"'workspace.filetype' must be a string, got {type(workspace['filetype']).__name__}.")
    else:
        workspace["filetype"] = workspace["filetype"].lower().strip()
        if workspace["filetype"] not in WORKSPACE_FILETYPES: error(f"'workspace.filetype' must be either 'tsv' or 'csv', got {workspace['filetype']!r}.")

    # ========== protein_markers ==========
    protein_markers = config.get("protein_markers")

    if protein_markers is None:
        error("'protein_markers' setting is missing.")
        protein_markers = []
    elif not isinstance(protein_markers, list):
        error(f"'protein_markers' must be a list, got {type(protein_markers).__name__}.")
        protein_markers = []
    else:
        if not protein_markers: error("At least one primary protein marker is required.")

        seen_markers: set[str] = set()

        for index, marker in enumerate(protein_markers):
            if not isinstance(marker, str):
                error(f"'protein_markers[{index}]' must be a string, got {type(marker).__name__}.")
                continue
            if not _valid_marker_name(marker): error(f"Protein marker {marker!r} is not valid.")
            if marker in seen_markers: error(f"Primary protein marker {marker!r} appears more than once in 'protein_markers'.")
            seen_markers.add(marker)

    primary_markers = {marker for marker in protein_markers if isinstance(marker, str)}

    # ========== cell_types ==========
    cell_types = config.get("cell_types")

    if cell_types is None:
        error("'cell_types' setting is missing.")
        cell_types = {}
    elif not isinstance(cell_types, dict):
        error(f"'cell_types' must be a dictionary, got {type(cell_types).__name__}.")
        cell_types = {}
    else:
        if not cell_types: error("At least one cell type is required.")

        for cell_type, definition in cell_types.items():
            if not isinstance(cell_type, str) or not cell_type.strip(): error(f"Cell type name {cell_type!r} must be a non-empty string.")

            if not isinstance(definition, list):
                error(f"Definition for cell type {cell_type!r} must be a list, got {type(definition).__name__}.")
                continue

            if len(definition) != len(protein_markers): error(f"Cell type {cell_type!r} contains {len(definition)} protein-marker values, but 'protein_markers' contains {len(protein_markers)} markers.")

            for index, value in enumerate(definition):
                if _is_binary_or_null(value): continue

                marker = (protein_markers[index] if index < len(protein_markers) else f"index {index}")
                error(f"Value {value!r} for protein marker {marker!r}, defining {cell_type!r}, must be 0, 1, or null.")

    # ========== sub_cell_types ==========
    if "sub_cell_types" not in config or config["sub_cell_types"] is None: config["sub_cell_types"] = {}

    sub_cell_types = config["sub_cell_types"]
    child_to_parents: dict[str, list[str]] = defaultdict(list)
    edge_markers: dict[tuple[str, str], set[str]] = {}
    all_children: set[str] = set()

    if not isinstance(sub_cell_types, dict):
        error(f"'sub_cell_types' must be a dictionary, got {type(sub_cell_types).__name__}.")
        sub_cell_types = {}
        config["sub_cell_types"] = {}
    else:
        for parent, children in sub_cell_types.items():
            if not isinstance(parent, str) or not parent.strip(): error(f"Subtype parent {parent!r} must be a non-empty string.")
            if not isinstance(children, dict):
                error(f"Descendants of {parent!r} must be a dictionary, got {type(children).__name__}.")
                continue

            if not children: error(f"Subtype parent {parent!r} must contain at least one descendant.")

            for child, definition in children.items():
                if isinstance(child, str):
                    all_children.add(child)
                    if isinstance(parent, str): child_to_parents[child].append(parent)

                if not isinstance(child, str) or not child.strip(): error(f"Cell subtype {child!r} must be a non-empty string.")

                if child == parent: error(f"Cell subtype {child!r} cannot be its own parent.")

                if not isinstance(definition, dict):
                    error(f"Protein-marker definition for {child!r} must be a dictionary, got {type(definition).__name__}.")
                    continue

                if not definition: error(f"At least one protein marker is required to define subtype {child!r}.")

                current_markers: set[str] = set()

                for marker, value in definition.items():
                    if not isinstance(marker, str):
                        error(f"Marker names defining {child!r} must be strings, got {type(marker).__name__}.")
                        continue

                    current_markers.add(marker)

                    if not _valid_marker_name(marker): error(f"Subtype marker {marker!r}, defining {child!r}, does not have a valid protein-marker name.")

                    if marker in primary_markers: error(f"Primary protein marker {marker!r} cannot also be used to define subtype {child!r}.")

                    # Subtype markers actively define the subtype, so they must be 0 or 1 rather than null.
                    if not _is_binary(value): error(f"Value {value!r} for subtype marker {marker!r}, defining {child!r}, must be 0 or 1.")

                if isinstance(parent, str) and isinstance(child, str): edge_markers[(parent, child)] = current_markers

        root_cell_types = set(cell_types)
        known_nodes = root_cell_types | all_children

        for parent in sub_cell_types:
            if parent not in known_nodes: error(f"Subtype parent {parent!r} is neither a primary cell type nor a defined cell subtype.")

        for child, parents in child_to_parents.items():
            if len(parents) > 1: error(f"Cell subtype {child!r} has multiple parents: {parents}.")

        graph = { parent:list(children) for parent, children in sub_cell_types.items() if isinstance(parent, str) and isinstance(children, dict)}

        # ---------- Detect cycles ----------
        node_state: dict[str, int] = {}
        cycle_found = False

        def detect_cycles(node: str, path: list[str]) -> None:
            nonlocal cycle_found
            node_state[node] = 1

            for child in graph.get(node, []):
                if not isinstance(child, str): continue
                if node_state.get(child) == 1:
                    cycle_found = True
                    error("Cycle detected in 'sub_cell_types': " + " -> ".join(path + [child]))
                elif node_state.get(child, 0) == 0: detect_cycles(child, path + [child])
            node_state[node] = 2

        for node in graph:
            if node_state.get(node, 0) == 0: detect_cycles(node, [node])

        # ---------- Detect orphan branches ----------
        reachable = set(root_cell_types)
        changed = True

        while changed:
            changed = False
            for parent, children in graph.items():
                if parent not in reachable: continue

                for child in children:
                    if isinstance(child, str) and child not in reachable:
                        reachable.add(child)
                        changed = True

        for parent in sub_cell_types:
            if parent not in reachable and parent in known_nodes: error(f"Subtype branch starting at {parent!r} cannot be reached from any primary cell type.")

        # ---------- Ancestor protein-marker rule ----------
        # Markers used by ancestors cannot be reused by descendants.
        # Sibling branches are validated independently, so siblings may use
        # the same marker when they share the same parent.
        if not cycle_found:

            def validate_branch(
                parent: str,
                ancestor_markers: set[str],
                path: list[str],
            ) -> None:
                for child in graph.get(parent, []):
                    if not isinstance(child, str): continue

                    markers = edge_markers.get((parent, child), set())
                    repeated = markers & ancestor_markers

                    if repeated: error(f"Subtype {child!r} reuses protein marker(s) already used by an ancestor: {sorted(repeated)}. Path: {' -> '.join(path + [child])}.")

                    validate_branch(child, ancestor_markers | markers, path + [child],)

            for root in root_cell_types: validate_branch(root, set(primary_markers), [root])

    # ========== spatial_plot ==========
    spatial_plot = config.get("spatial_plot")

    if spatial_plot is None: error("'spatial_plot' setting is missing.")
    elif not isinstance(spatial_plot, dict): error(f"'spatial_plot' must be a dictionary, got {type(spatial_plot).__name__}.")
    else:
        # The supplied YAML uses "plot_options". "spatial_options" is
        # accepted as a backwards-compatible alias.
        if "plot_options" not in spatial_plot and "spatial_options" in spatial_plot:
            spatial_plot["plot_options"] = spatial_plot.pop("spatial_options")
            warning("'spatial_options' was renamed to 'plot_options' to match the configuration schema.")

        if "plot_options" not in spatial_plot:
            warning("'spatial_plot.plot_options' is missing. Defaulting to 'all'.")
            spatial_plot["plot_options"] = "all"

        plot_option = spatial_plot["plot_options"]

        if not isinstance(plot_option, str) or plot_option not in PLOT_OPTIONS: error(f"'spatial_plot.plot_options' must be one of {sorted(PLOT_OPTIONS)}, got {plot_option!r}.")

        if "custom_plot" not in spatial_plot or spatial_plot["custom_plot"] is None: spatial_plot["custom_plot"] = []

        custom_plot = spatial_plot["custom_plot"]

        if not isinstance(custom_plot, list): error(f"'spatial_plot.custom_plot' must be a list, got {type(custom_plot).__name__}.")
        else:
            if plot_option == "custom" and not custom_plot:
                warning("'spatial_plot.plot_options' is 'custom', but 'custom_plot' is empty. Changing plot_options to 'none'.")
                spatial_plot["plot_options"] = "none"
            elif plot_option != "custom" and custom_plot: warning("'spatial_plot.custom_plot' contains entries, but plot_options is not 'custom'. Custom entries will not normally be used.")

            seen_files: set[str] = set()
            filetype = workspace.get("filetype")

            for index, item in enumerate(custom_plot):
                if not isinstance(item, dict):
                    error(f"'spatial_plot.custom_plot[{index}]' must be a dictionary, got {type(item).__name__}.")
                    continue

                expected_keys = {"fileName", "type"}
                if set(item) != expected_keys: error(f"'spatial_plot.custom_plot[{index}]' must contain exactly the keys 'fileName' and 'type'.")

                file_name = item.get("fileName")
                custom_type = item.get("type")

                if not isinstance(file_name, str) or not file_name.strip(): error(f"'spatial_plot.custom_plot[{index}].fileName' must be a non-empty string.")
                else:
                    file_name = file_name.strip()
                    if "_objects" not in Path(file_name).stem: error(f"Custom plot file {file_name!r} must contain '_objects' in its filename.")
                    if filetype in WORKSPACE_FILETYPES:
                        expected_extension = f".{filetype}"
                        if Path(file_name).suffix.lower() != expected_extension: error(f"Custom plot file {file_name!r} must use {expected_extension!r} to match 'workspace.filetype'.")
                        else:
                            if input_dir_checked is not None and not Path(input_dir_checked / file_name).is_file(): error(f"Custom plot file {file_name!r} does not exist or is misspelled.")
                    if file_name in seen_files: error(f"Custom plot file {file_name!r} appears more than once.")
                    seen_files.add(file_name)

                if (not isinstance(custom_type, str) or custom_type not in CUSTOM_PLOT_TYPES):
                    error(f"'spatial_plot.custom_plot[{index}].type' must be one of {sorted(CUSTOM_PLOT_TYPES)}, got {custom_type!r}.")

        if "dpi" in spatial_plot:
            dpi = spatial_plot["dpi"]
            if (isinstance(dpi, bool) or not isinstance(dpi, (int, float)) or dpi <= 0): error("'spatial_plot.dpi' must be a positive number.")

        if "scatter_point_size" in spatial_plot:
            point_size = spatial_plot["scatter_point_size"]
            if (isinstance(point_size, bool) or not isinstance(point_size, (int, float)) or point_size <= 0): error("'spatial_plot.scatter_point_size' must be a positive number.")

    # ========== custom_colors ==========
    custom_colors = config.get("custom_colors")
    hierarchy_names = set(cell_types) | all_children | set(['Other cells'])

    if custom_colors is None: error("'custom_colors' setting is missing.")
    elif not isinstance(custom_colors, dict): error(f"'custom_colors' must be a dictionary, got {type(custom_colors).__name__}.")
    else:
        # if 'Other cells' not in custom_colors.keys(): error(f"'Other cells' is missing in 'custom_colors' keys.")
        for cell_name, color in custom_colors.items():
            if not isinstance(cell_name, str) or not cell_name.strip(): error("Keys in 'custom_colors' must be non-empty strings.")

            if not isinstance(color, str) or not HEX_COLOR_RE.fullmatch(color): error(f"Color for {cell_name!r} must be a hexadecimal color in the format '#RRGGBB', got {color!r}.")

        missing_colors = hierarchy_names - set(custom_colors)
        if missing_colors: error("Missing custom color(s) for: " + ", ".join(sorted(missing_colors)) + ".")

        extra_colors = set(custom_colors) - hierarchy_names
        if extra_colors: warning("'custom_colors' contains name(s) that are not present in the cell hierarchy: " + ", ".join(sorted(extra_colors)) + ".")

        color_to_cells: dict[str, list[str]] = defaultdict(list)
        for cell_name, color in custom_colors.items():
            if isinstance(color, str) and HEX_COLOR_RE.fullmatch(color): color_to_cells[color.lower()].append(cell_name)

        for color, cells in color_to_cells.items():
            if len(cells) > 1: warning(f"Color {color!r} is used by multiple cell types/subtypes: {cells}.")

    # ========== Boolean settings present in the supplied configuration ==========
    for key in ("locally_save_anndata_files", "overwrite_existing_files"):
        if key in config and type(config[key]) is not bool: error(f"'{key}' must be true or false, got {config[key]!r}.")

    # ========== Defaults + arrows_orientation validation ==========
    if "arrows_orientation" not in config: config["arrows_orientation"] = "contains"

    arrows_orientation = config["arrows_orientation"]

    if not isinstance(arrows_orientation, str): 
        warning(f"'arrows_orientation' must be a string, got {type(arrows_orientation).__name__}. Therefore, 'arrows_orientation' is set to default value.")
        config["arrows_orientation"] = "contains"
    elif arrows_orientation not in ARROWS_ORIENTATIONS: 
        warning(f"'arrows_orientation' must be one of {sorted(ARROWS_ORIENTATIONS)}, got {arrows_orientation!r}. Therefore, 'arrows_orientation' is set to default value.")
        config["arrows_orientation"] = "contains"

    # ========== Final result ==========
    if errors: raise ConfigValidationError(errors)
    return config


def _resolve_workspace_path(path_value: str, config_path: Path) -> Path:
    """
    Resolve a configured workspace path to an absolute path.

    Parameters
    ----------
    path_value : str
        User-provided input or output directory.
    config_path : pathlib.Path
        Path to the YAML configuration file.

    Returns
    -------
    pathlib.Path
        Absolute normalized workspace path.

    Notes
    -----
    Relative workspace paths are interpreted relative to the directory
    containing the YAML file rather than the current working directory.
    """
    path = Path(path_value).expanduser()
    if not path.is_absolute(): path = config_path.parent / path
    return path.resolve()


def load_config(config_path: str | Path) -> dict:
    """Load, validate, and normalize a SPPRAN YAML configuration file.

    The function performs four steps:

    1. Validate the configuration-file path.
    2. Parse YAML safely while rejecting duplicate keys.
    3. Validate and normalize the configuration schema.
    4. Convert workspace paths to ``Path`` objects and verify that the
       user-defined input and output directories exist.

    Relative ``input_dir`` and ``output_dir`` values are resolved relative
    to the YAML file's directory.

    Parameters
    ----------
    config_path : str | Path
        Path to the YAML configuration file.

    Returns
    -------
    dict
        Validated and normalized configuration dictionary. The returned
        ``workspace.input_dir`` and ``workspace.output_dir`` values are
        absolute ``Path`` objects.

    Raises
    ------
    FileNotFoundError
        If the configuration file, input directory, or output directory
        does not exist.
    IsADirectoryError
        If ``config_path`` points to a directory instead of a file.
    NotADirectoryError
        If ``input_dir`` or ``output_dir`` exists but is not a directory.
    ConfigValidationError
        If the YAML is syntactically valid but contains invalid settings.
    ValueError
        If the YAML file cannot be parsed.
    """
    config_path = Path(config_path).expanduser()
    if not config_path.exists(): raise FileNotFoundError(f"Configuration file not found: '{config_path}'.")
    if not config_path.is_file(): raise IsADirectoryError(f"Configuration path '{config_path}' exists, but it is not a file.")

    try:
        with config_path.open("r", encoding="utf-8") as file:
            config = yaml.load(file, Loader=UniqueKeyLoader)
    except yaml.YAMLError as exc:
        raise ValueError(f"There's an error parsing the configuration file '{config_path}'. {REVIEW_YAML_MESSAGE}") from exc

    config = _check_yaml(config)

    input_dir = _resolve_workspace_path(config["workspace"]["input_dir"], config_path,)
    output_dir = _resolve_workspace_path(config["workspace"]["output_dir"], config_path,)
    config["workspace"]["input_dir"] = input_dir
    config["workspace"]["output_dir"] = output_dir
    # if not input_dir.exists(): raise FileNotFoundError(f"Input directory '{input_dir}' does not exist or is misspelled.")
    # if not input_dir.is_dir(): raise NotADirectoryError(f"Input path '{input_dir}' exists, but it is not a directory.")

    # # The output directory is intentionally user-defined and is not created automatically. This keeps result organization explicit and predictable.
    # if not output_dir.exists(): raise FileNotFoundError(f"Output directory '{output_dir}' does not exist or is misspelled. Please create the desired output directory and specify it in the configuration file.")
    # if not output_dir.is_dir(): raise NotADirectoryError(f"Output path '{output_dir}' exists, but it is not a directory.")
    return config

__all__ = ["ConfigValidationError", "load_config",]
