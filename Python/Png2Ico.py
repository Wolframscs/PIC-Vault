import pandas as pd
import numpy as np
import datetime
import tkinter as tk
from tkinter import filedialog, messagebox

# 主窗口
def main():
    def run():
        file = file_entry.get()
        sheet = sheet_entry.get()
        if not file:
            messagebox.showerror('Error', 'Please select an Excel file!')
            return
        if not sheet:
            messagebox.showerror('Error', 'Please enter the sheet name!')
            return
        try:
            df_input = pd.read_excel(file, sheet_name=sheet)
            if 'SOC' not in df_input.columns or 'Rate' not in df_input.columns:
                msg = f"Sheet '{sheet}' in Excel file '{file}' must contain both 'SOC' and 'Rate' columns!"
                msg += f" Actual columns: {list(df_input.columns)}"
                raise ValueError(msg)
            soc = df_input['SOC'].tolist()
            rate = df_input['Rate'].tolist()
            time_list = []
            rate_list = []
            soc_list = []
            current_time = 0
            current_soc = 0.0
            time_list.append(current_time)
            soc_list.append(current_soc)
            rate_list.append(rate[1])
            while current_soc < 1:
                for idx in range(len(soc) - 1):
                    if soc[idx] <= current_soc < soc[idx + 1]:
                        r = rate[idx + 1]
                        break
                    elif current_soc >= soc[-1]:
                        r = rate[-1]
                        break
                current_time += 1
                current_soc += r * (1 / 3600)
                if current_soc > 1:
                    current_soc = 1
                time_list.append(current_time)
                soc_list.append(current_soc)
                rate_list.append(r)
            df = pd.DataFrame({'Time': time_list, 'SOC': soc_list, 'Rate': rate_list})
            now = datetime.datetime.now()
            filename = now.strftime('Time_Rate_%Y%m%d_%H%M%S.xlsx')
            df.to_excel(filename, index=False)
            messagebox.showinfo('Success', f'File saved as {filename}')
        except Exception as e:
            messagebox.showerror('Error', str(e))

    def open_file():
        file_path = filedialog.askopenfilename(filetypes=[('Excel Files', '*.xlsx')])
        if file_path:
            file_entry.delete(0, tk.END)
            file_entry.insert(0, file_path)

    root = tk.Tk()
    root.title('SOC Rate EXE')
    root.geometry('400x180')

    tk.Label(root, text='Excel file:').grid(row=0, column=0, padx=10, pady=10, sticky='e')
    file_entry = tk.Entry(root, width=30)
    file_entry.grid(row=0, column=1, padx=5)
    tk.Button(root, text='Open', command=open_file).grid(row=0, column=2, padx=5)

    tk.Label(root, text='Sheet name:').grid(row=1, column=0, padx=10, pady=10, sticky='e')
    sheet_entry = tk.Entry(root, width=30)
    sheet_entry.insert(0, 'Sheet1')
    sheet_entry.grid(row=1, column=1, padx=5)

    tk.Button(root, text='Run', command=run, width=15).grid(row=2, column=1, pady=20)

    root.mainloop()

if __name__ == '__main__':
    main()
