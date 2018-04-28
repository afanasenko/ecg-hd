#!/usr/bin/env bash

# locate ui compiler
UIC=/anaconda2/bin/pyuic5

if [ -f $UIC ]
then
    $UIC heartapp_mainwindow.ui > ui_mainwindow.py
else
    echo "pyuic not found"
fi

