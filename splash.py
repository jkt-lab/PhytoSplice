import tkinter as tk
from tkinter import ttk
import sys
import os

def main():
    try:
        root = tk.Tk()
        root.overrideredirect(True) # Frameless
        
        # Center the window
        width = 400
        height = 180
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        x = (screen_width // 2) - (width // 2)
        y = (screen_height // 2) - (height // 2)
        root.geometry(f'{width}x{height}+{int(x)}+{int(y)}')
        
        # Colors matching the dashboard
        bg_color = '#222d32'
        text_color = '#ffffff'
        accent_color = '#3c8dbc'
        
        root.configure(bg=bg_color)
        
        # Main Frame
        frame = tk.Frame(root, bg=bg_color)
        frame.pack(expand=True, fill='both', padx=20, pady=20)
        
        # App Title
        label = tk.Label(frame, text="PhytoSplice", font=("Helvetica", 16, "bold"), fg=text_color, bg=bg_color)
        label.pack(pady=(10, 5))
        
        # Subtitle
        status = tk.Label(frame, text="Starting application server...", font=("Helvetica", 10), fg="#b8c7ce", bg=bg_color)
        status.pack(pady=(0, 20))
        
        # Progress bar
        style = ttk.Style()
        style.theme_use('default')
        style.configure("Custom.Horizontal.TProgressbar", background=accent_color, troughcolor="#1a2226", bordercolor=bg_color, thickness=6)
        
        pb = ttk.Progressbar(frame, style="Custom.Horizontal.TProgressbar", orient="horizontal", length=350, mode="indeterminate")
        pb.pack(pady=10)
        pb.start(15)
        
        root.mainloop()
    except Exception as e:
        # Fallback if GUI fails (e.g. no display), just exit silent
        print(f"Splash error: {e}")
        sys.exit(0)

if __name__ == "__main__":
    main()
