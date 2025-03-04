import aeolis
from tkinter import *
from tkinter import ttk

configfile = r'C:\Users\svries\Documents\GitHub\OE_aeolis-python\aeolis\examples\2D\Barchan_dune\aeolis.txt'

dic = aeolis.inout.read_configfile(configfile)

# Create a window
root = Tk()
root.title('Aeolis')
root.geometry('400x400')

# Create a frame
frame = ttk.Frame(root)
frame.pack()

# create tabs
tab_control = ttk.Notebook(root)
tab1 = ttk.Frame(tab_control)
tab2 = ttk.Frame(tab_control)
tab3 = ttk.Frame(tab_control)

# Helper function to create label and entry
def create_label_entry(tab, text, value, row):
    label = ttk.Label(tab, text=text)
    label.grid(row=row, column=0, sticky=W)
    entry = ttk.Entry(tab)
    entry.insert(0, str(value))
    entry.grid(row=row, column=1, sticky=W)
    return entry

## this is tab1
# Display in tab1 xgrid_file, ygrid_file, bed_file, ne_file, veg_file keys with their associated values
xgrid_file = dic.get('xgrid_file', '')
ygrid_file = dic.get('ygrid_file', '')
bed_file = dic.get('bed_file', '')
ne_file = dic.get('ne_file', '')
veg_file = dic.get('veg_file', '')
threshold_file = dic.get('threshold_file', '')
fence_file = dic.get('fence_file', '')
wave_mask = dic.get('wave_mask', '')
tide_mask = dic.get('tide_mask', '')
threshold_mask = dic.get('threshold_mask', '')

xgrid_file_entry = create_label_entry(tab1, "xgrid_file:", xgrid_file, 0)
ygrid_file_entry = create_label_entry(tab1, "ygrid_file:", ygrid_file, 1)
bed_file_entry = create_label_entry(tab1, "bed_file:", bed_file, 2)
ne_file_entry = create_label_entry(tab1, "ne_file:", ne_file, 3)
veg_file_entry = create_label_entry(tab1, "veg_file:", veg_file, 4)
threshold_file_entry = create_label_entry(tab1, "threshold_file:", threshold_file, 5)
fence_file_entry = create_label_entry(tab1, "fence_file:", fence_file, 6)
wave_mask_entry = create_label_entry(tab1, "wave_mask:", wave_mask, 7)
tide_mask_entry = create_label_entry(tab1, "tide_mask:", tide_mask, 8)
threshold_mask_entry = create_label_entry(tab1, "threshold_mask:", threshold_mask, 9)

tab_control.add(tab1, text='Domain')

# Lets make tab 2
# Display in tab2 the keys: tstart, tstop, dt, restart, refdate
tstart = dic.get('tstart', '')
tstop = dic.get('tstop', '')
dt = dic.get('dt', '')
restart = dic.get('restart', '')
refdate = dic.get('refdate', '')

tstart_entry = create_label_entry(tab2, "tstart:", tstart, 0)
tstop_entry = create_label_entry(tab2, "tstop:", tstop, 1)
dt_entry = create_label_entry(tab2, "dt:", dt, 2)
restart_entry = create_label_entry(tab2, "restart:", restart, 3)
refdate_entry = create_label_entry(tab2, "refdate:", refdate, 4)

tab_control.add(tab2, text='Timeframe')

# lets make tab 3

tab_control.add(tab3, text='Sediment Transport')
tab_control.pack(expand=1, fill='both')


# # Create a button
# button = ttk.Button(frame, text='Run')
# button.pack()

# # Create a text box
# text = Text(frame, width=40, height=10)
# text.pack()

# add a button that saves the updated values in a text file
def save(*args):
    dic['xgrid_file'] = xgrid_file_entry.get()
    dic['ygrid_file'] = ygrid_file_entry.get()
    dic['bed_file'] = bed_file_entry.get()
    dic['ne_file'] = ne_file_entry.get()
    dic['veg_file'] = veg_file_entry.get()
    dic['threshold_file'] = threshold_file_entry.get()
    dic['fence_file'] = fence_file_entry.get()
    dic['wave_mask'] = wave_mask_entry.get()
    dic['tide_mask'] = tide_mask_entry.get()
    dic['threshold_mask'] = threshold_mask_entry.get()

    dic['tstart'] = tstart_entry.get()
    dic['tstop'] = tstop_entry.get()
    dic['dt'] = dt_entry.get()
    dic['restart'] = restart_entry.get()
    dic['refdate'] = refdate_entry.get()

    aeolis.inout.write_configfile(configfile + '2', dic)

    print('Saved!')

save_button = ttk.Button(tab3, text='Save', command=save)
save_button.pack()





# display gui
root.mainloop()




# print(dic)