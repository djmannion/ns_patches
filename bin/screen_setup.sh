#!/usr/bin/env bash

xrandr --output VGA2 --primary
./dummy_window
xrandr --output eDP2 --off
./dummy_window
xrandr --output VGA2 --primary
./dummy_window
xrandr --output eDP2 --mode 1024x768
./dummy_window
