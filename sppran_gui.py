import tkinter as tk
from tkinter import Tk, Button, Label, filedialog, messagebox, scrolledtext
from spatial_proteomics.core import main_gui

def choose_config_file():
    """Opens a file dialog to select a YAML config file and displays its content in the text area."""
    config_file_path = filedialog.askopenfilename(
        title="Select config YAML file",
        filetypes=[("YAML files", "*.yaml, *.yml")]
    )
    config_label.config(text=config_file_path)
    with open(config_file_path, "r", encoding="utf-8", errors="replace") as file:
        content = file.read()
    text_area.delete(1.0, tk.END)
    text_area.insert(tk.END, content)
    root.config_path = config_file_path

def run_pipeline():
    """Runs the SpPrAn analysis pipeline using the selected config file."""
    if not getattr(root, 'config_path', None):
        messagebox.showerror("Missing config file", "Please select a config YAML file.")
        return
    try:
        main_gui(root.config_path)
        messagebox.showinfo("Done", "SpPrAn analysis pipeline finished.")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

root = Tk()
root.title("SpPrAn")

Button(root, text="Choose config YAML file", command=choose_config_file).pack(pady=10)
config_label = Label(root, text="No config file selected")
config_label.pack()

text_area = scrolledtext.ScrolledText(root, wrap=tk.WORD, font=("Courier", 10))
text_area.pack(expand=True, fill="both", padx=5, pady=5)

Button(root, text="Run analysis", command=run_pipeline).pack(pady=10)
root.mainloop()
