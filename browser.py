#!/usr/bin/env python
# coding: utf-8

from Tkinter import *
import tkFileDialog
import os

import crossplot


class App:
    def __init__(self, master):
        self.master = master

        # флаги для расчета графиков
        self.show_scope = BooleanVar(value=False)

        self.show_freq = BooleanVar(value=True)

        self.show_hist = BooleanVar(value=True)

        self.show_cwt = BooleanVar(value=True)

        self.fixbase = BooleanVar(value=True)
        self.mains_filter = BooleanVar(value=False)

        self.offset_1 = IntVar(value=0)
        self.offset_2 = IntVar(value=0)
        self.processing_window = IntVar(value=60)

        self.writefile_base = StringVar(value=u"опыт_0001")

        # call start to initialize to create the UI elemets
        self.start()

    def start(self):
        self.master.title(u"Выбор файлов для анализа")

        # CREATE A TEXT/LABEL
        # then, put in the first row (row=0) and in the 2nd column (column=1), align it to "West"/"W"
        Label(self.master, text="1:").grid(row=0, column=0)

        # CREATE A TEXTBOX
        self.filename1 = Entry(self.master, width=64)
        self.filename1.focus_set()
        self.filename1.grid(row=0, column=1)

        # CREATE A BUTTON WITH "ASK TO OPEN A FILE"
        open_file_a = Button(self.master, text=u"Выбор файла А",
                                command=self.browse_file_1)
        open_file_a.grid(row=0, column=2)

        Label(self.master, text=u"2:").grid(row=1, column=0)

        # CREATE A TEXTBOX
        self.filename2 = Entry(self.master, width=64)
        self.filename2.grid(row=1, column=1)

        # CREATE A BUTTON WITH "ASK TO OPEN A FILE"
        open_file_b = Button(self.master, text=u"Выбор файла Б",
                                command=self.browse_file_2)
        open_file_b.grid(row=1, column=2)

        chk1 = Checkbutton(self.master, text=u"Осциллограммы",
                           var=self.show_scope)
        chk1.grid(row=2, column=1, sticky=W)

        chk2 = Checkbutton(self.master, text=u"Спектры",
                           var=self.show_freq)
        chk2.grid(row=3, column=1, sticky=W)

        chk3 = Checkbutton(self.master, text=u"Гистограммы",
                           var=self.show_hist)
        chk3.grid(row=4, column=1, sticky=W)

        chk4 = Checkbutton(self.master, text=u"Вейвлеты",
                           var=self.show_cwt)
        chk4.grid(row=5, column=1, sticky=W)

        Label(self.master, text=u"Параметры обработки:").grid(row=6, column=1)

        chk5 = Checkbutton(self.master, text=u"Подавление сетевой помехи",
                           var=self.mains_filter)
        chk5.grid(row=7, column=1, sticky=W)

        chk6 = Checkbutton(self.master, text=u"Подавление дрейфа изолинии",
                           var=self.fixbase)

        chk6.grid(row=8, column=1, sticky=W)

        Label(self.master, text=u"Начало обработки в файле А, с:").grid(
            row=9, column=1, sticky=W)
        Entry(self.master, textvar=self.offset_1).grid(
            row=9, column=1, sticky=E)

        Label(self.master, text=u"Начало обработки в файле Б, с:").grid(
            row=10, column=1, sticky=W)
        Entry(self.master, textvar=self.offset_2).grid(
            row=10, column=1, sticky=E)

        Label(self.master, text=u"Обрабатываемый интервал, с:").grid(
            row=11, column=1, sticky=W)
        Entry(self.master, textvar=self.processing_window).grid(
            row=11, column=1, sticky=E)

        # CREATE RADIO BUTTONS
        RADIO_BUTTON = [
            (W, u"Сохранить графики на диске", "disk"),
            (E, u"Показать графики на экране", "screen")
        ]

        Label(self.master, text=u"Вывод результата:").grid(row=12, column=1)

        self.radio_var = StringVar()
        self.radio_var.set("screen")

        for i, text, item in RADIO_BUTTON:
            # setup each radio button. variable is set to the self.radio_var
            # and the value is set to the "item" in the for loop
            self.radio = Radiobutton(self.master, text=text,
                                     variable=self.radio_var, value=item)
            self.radio.grid(row=13, column=1, sticky=i)

        Entry(self.master, textvar=self.writefile_base).grid(
            row=14, column=1, sticky=W)

        # now for a button
        self.submit = Button(self.master, text="Рассчитать",
                             command=self.start_processing, width=12, height=3)
        self.submit.grid(row=15, column=2)

    def start_processing(self):
        crossplot.processing(
            file1=self.filename1.get(),
            file2=self.filename2.get(),
            use_mains_filter=self.mains_filter.get(),
            show_scope=self.show_scope.get(),
            show_freq=self.show_freq.get(),
            show_hist=self.show_hist.get(),
            show_cwt=self.show_cwt.get(),
            offset_first=self.offset_1.get(),
            offset_second=self.offset_2.get(),
            window=self.processing_window.get(),
            file_output=self.writefile_base.get() if self.radio_var.get() ==
                                                     "disk" else ""
        )
        print(u"Готово!")

    def browse_file_1(self):

        f = self.filename1.get()
        if not f:
            f = __file__

        filename = tkFileDialog.askopenfilename(
            title=u"Выбор первого файла",
            filetypes=(("Zetlab", "*.ana"), ("all files", "*.*")),
            initialdir=os.path.dirname(f)
        )
        self.filename1.delete(0, END)
        self.filename1.insert(0, filename)

    def browse_file_2(self):
        f = self.filename2.get()
        if not f:
            f = __file__

        filename = tkFileDialog.askopenfilename(
            title=u"Выбор второго файла",
            filetypes=(("Zetlab", "*.ana"), ("all files", "*.*")),
            initialdir=os.path.dirname(f)
        )
        self.filename2.delete(0, END)
        self.filename2.insert(0, filename)


root = Tk()
app = App(root)
root.mainloop()