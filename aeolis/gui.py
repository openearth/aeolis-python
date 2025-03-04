import aeolis
from tkinter import *
from tkinter import ttk, filedialog
import os

configfile = r'C:\Users\svries\Documents\GitHub\OE_aeolis-python\aeolis\examples\2D\Barchan_dune\aeolis.txt'

def prompt_file():
    file_path = filedialog.askopenfilename(
        initialdir=os.path.dirname(configfile),
        title="Select config file",
        filetypes=(("Text files", "*.txt"), ("All files", "*.*"))
    )
    return file_path if file_path else configfile

configfile = prompt_file()
dic = aeolis.inout.read_configfile(configfile)

class AeolisGUI:
    def __init__(self, root, dic):
        self.root = root
        self.dic = dic
        self.root.title('Aeolis')
        self.create_widgets()

    def create_widgets(self):
        tab_control = ttk.Notebook(self.root)
        self.create_domain_tab(tab_control)
        self.create_timeframe_tab(tab_control)
        self.create_boundary_conditions_tab(tab_control)
        self.create_sediment_transport_tab(tab_control)
        tab_control.pack(expand=1, fill='both')

    def create_label_entry(self, tab, text, value, row):
        label = ttk.Label(tab, text=text)
        label.grid(row=row, column=0, sticky=W)
        entry = ttk.Entry(tab)
        entry.insert(0, str(value))
        entry.grid(row=row, column=1, sticky=W)
        return entry

    def create_domain_tab(self, tab_control):
        tab1 = ttk.Frame(tab_control)
        tab_control.add(tab1, text='Domain')

        fields = ['xgrid_file', 'ygrid_file', 'bed_file', 'ne_file', 'veg_file', 'threshold_file', 'fence_file', 'wave_mask', 'tide_mask', 'threshold_mask']
        self.entries = {field: self.create_label_entry(tab1, f"{field}:", self.dic.get(field, ''), i) for i, field in enumerate(fields)}

        fig_frame = ttk.Frame(tab1)
        fig_frame.grid(row=10, column=0, columnspan=2, pady=10)
        fig_label = ttk.Label(fig_frame, text="Figures:")
        fig_label.pack()

        fig_canvas_frame = ttk.Frame(tab1)
        fig_canvas_frame.grid(row=0, column=2, rowspan=10, padx=10, pady=10, sticky=N)
        self.fig_canvas = Canvas(fig_canvas_frame, width=300, height=200, bg='white')
        self.fig_canvas.pack()

        update_button = ttk.Button(fig_canvas_frame, text="Update Figure", command=self.update_figure)
        update_button.pack()

    def create_timeframe_tab(self, tab_control):
        tab2 = ttk.Frame(tab_control)
        tab_control.add(tab2, text='Timeframe')

        fields = ['tstart', 'tstop', 'dt', 'restart', 'refdate']
        self.entries.update({field: self.create_label_entry(tab2, f"{field}:", self.dic.get(field, ''), i) for i, field in enumerate(fields)})

    def create_boundary_conditions_tab(self, tab_control):
        tab3 = ttk.Frame(tab_control)
        tab_control.add(tab3, text='Boundary Conditions')

        fields = ['boundary1', 'boundary2', 'boundary3']
        self.entries.update({field: self.create_label_entry(tab3, f"{field}:", self.dic.get(field, ''), i) for i, field in enumerate(fields)})

    def create_sediment_transport_tab(self, tab_control):
        tab4 = ttk.Frame(tab_control)
        tab_control.add(tab4, text='Sediment Transport')

        save_button = ttk.Button(tab4, text='Save', command=self.save)
        save_button.pack()

    def update_figure(self):
        self.fig_canvas.create_rectangle(50, 50, 250, 150, fill="blue")

    def save(self):
        for field, entry in self.entries.items():
            self.dic[field] = entry.get()
        aeolis.inout.write_configfile(configfile + '2', self.dic)
        print('Saved!')

if __name__ == "__main__":
    root = Tk()
    app = AeolisGUI(root, dic)
    root.mainloop()
