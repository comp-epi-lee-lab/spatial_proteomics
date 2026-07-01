# Version 0.2
from email.mime import message
import os
import sys
from pathlib import Path

def get_app_cache_dir(app_name="SpPrAn"):
    """Returns the path to the application cache directory."""
    if sys.platform.startswith("win"):
        return Path(os.environ.get("LOCALAPPDATA", str(Path.home() / "AppData" / "Local"))) / app_name / "matplotlib"
    elif sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / app_name / "matplotlib"
    else:
        return Path.home() / ".cache" / app_name / "matplotlib"

mpl_cache_dir = get_app_cache_dir()
mpl_cache_dir.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(mpl_cache_dir)
os.environ["MPLBACKEND"] = "Agg"

import contextlib
import queue
import threading
import traceback
import tkinter as tk
from tkinter import Tk, Button, Label, filedialog, messagebox, scrolledtext
from spatial_proteomics.core import main_gui

class QueueWriter:
    """A file-like object that writes to a queue."""
    def __init__(self, queue):
        self.queue = queue

    def write(self, text):
        if text: self.queue.put(text)

    def flush(self):
        pass

class SpPrAnGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("SpPrAn")
        self.root.geometry("800x600")
        
        self.config_path = None
        self.queue = queue.Queue()

        self.main_frame = tk.Frame(root)
        self.log_frame = tk.Frame(root)
        
        self.build_main_frame()
        self.build_log_frame()

        self.show_main_frame()
        self.root.after(100, self.process_log_queue)

    def build_main_frame(self):
        Label(self.main_frame, text="SpPrAn: Spatial Proteomics Analysis", font=("Helvetica", 18, "bold")).pack(pady=20)        
        Button(self.main_frame, text="Choose config YAML file", command=self.choose_config_file, width=25).pack(pady=10)

        self.config_label = Label(self.main_frame, text="No config file selected", wraplength=600)
        self.config_label.pack(pady=10)
        self.text_area = scrolledtext.ScrolledText(self.main_frame, wrap=tk.WORD, height=22, width=85, font=("Courier", 10))
        self.text_area.pack(padx=5, pady=5, expand=True, fill=tk.BOTH)
        self.run_button = Button(self.main_frame, text="Run analysis", command=self.run_pipeline, width=25, state=tk.DISABLED)
        self.run_button.pack(pady=20)

    def build_log_frame(self):
        Label(self.log_frame, text="Run SpPrAn", font=("Helvetica", 16, "bold")).pack(pady=10)
        self.log_text = scrolledtext.ScrolledText(self.log_frame, wrap=tk.WORD, height=22, width=85, font=("Courier", 10))
        self.log_text.pack(padx=5, pady=5, expand=True, fill=tk.BOTH)
        self.done_button = Button(self.log_frame, text="Back to Main Menu", command=self.show_main_frame, state=tk.DISABLED, width=20)
        self.done_button.pack(pady=10)

    def show_main_frame(self):
        self.log_frame.pack_forget()
        self.main_frame.pack(expand=True, fill=tk.BOTH)

    def show_log_frame(self):
        self.main_frame.pack_forget()
        self.log_frame.pack(expand=True, fill=tk.BOTH)

    def choose_config_file(self):
        """Opens a file dialog to select a YAML config file and displays its content in the text area."""
        config_file_path = filedialog.askopenfilename(
            title="Select config YAML file",
            filetypes=[("YAML files", "*.yaml, *.yml")]
        )
        self.config_label.config(text=config_file_path)
        with open(config_file_path, "r", encoding="utf-8", errors="replace") as file:
            content = file.read()
        self.text_area.delete(1.0, tk.END)
        self.text_area.insert(tk.END, content)
        self.root.config_path = config_file_path
        self.run_button.config(state=tk.NORMAL)

    def run_pipeline(self):
        """Runs the SpPrAn analysis pipeline using the selected config file."""
        if not getattr(self.root, 'config_path', None):
            messagebox.showerror("Missing config file", "Please select a config YAML file.")
            return
        self.log_text.delete(1.0, tk.END)
        self.done_button.config(text="Running...", state=tk.DISABLED)
        self.show_log_frame()

        thread = threading.Thread(target=self.run_pipeline_thread, daemon=True)
        thread.start()

    def run_pipeline_thread(self):
        writer = QueueWriter(self.queue)
        try:
            with contextlib.redirect_stdout(writer), contextlib.redirect_stderr(writer):
                print("Starting SpPrAn analysis pipeline...")
                print(f"Using config file: {self.root.config_path}")
                main_gui(self.root.config_path)
                self.queue.put("SpPrAn analysis pipeline finished successfully.\n")
            self.queue.put("__DONE__")
        except Exception as e:
            error_message = f"[ERROR] An error occurred: {e}\n{traceback.format_exc()}"
            self.queue.put(error_message)
        
    def process_log_queue(self):
        try:
            while True:
                item = self.queue.get_nowait()
                if item == "__DONE__":
                    self.done_button.config(text="Back to Main Menu", state=tk.NORMAL)
                    messagebox.showinfo("Sp", "SpPrAn analysis pipeline finished.")
                    continue
                self.log_text.insert(tk.END, item)
                self.log_text.see(tk.END)
        except queue.Empty:
            pass
        self.root.after(100, self.process_log_queue)

if __name__ == "__main__":
    root = Tk()
    app = SpPrAnGUI(root)
    root.mainloop()
