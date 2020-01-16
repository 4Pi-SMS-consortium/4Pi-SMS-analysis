function [mephi]=coords2mephi(coords,msz)

mephi=coords./msz.*2.*pi-pi;