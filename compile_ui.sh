#!/usr/bin/env bash
# coding: utf-8
# Вызов компилятора графического интерфейса

# locate ui compiler
if [ $# -ge 1 ]
then
    UIC=$1
else
    UIC=/anaconda2/bin/pyuic5
fi

if [ -f $UIC ]
then
    $UIC heartapp_mainwindow.ui > ui_mainwindow.py
else
    echo "pyuic not found:"
    echo $UIC
fi


