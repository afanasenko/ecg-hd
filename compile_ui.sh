#!/usr/bin/env bash

# locate ui compiler
UIC=~/anaconda3/bin/pyuic4

$UIC heartapp_mainwindow.ui > ui_mainwindow.py
