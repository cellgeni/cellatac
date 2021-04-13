#!/usr/bin/python3

import os, sys
import argparse
import csv
import re

my_parser = argparse.ArgumentParser()

my_parser.add_argument("-i", "--input", help="the file to be searched through (default STDIN)")
my_parser.add_argument("-s", "--input_field_sep", default='\t', help="the character the file columns are separated by")
my_parser.add_argument("-c", "--colname", help="the name of the column to check value of")
my_parser.add_argument("-f", "--field_value", help="the value to use to get the data")
my_parser.add_argument("-o", "--output_field_sep", default='\t', help="the output file delimiter")
my_parser.add_argument("-H", action="store_false", help="should header be printed")

args = my_parser.parse_args()

if args.colname is None:
  print('Error: Column name is not specified', file=sys.stderr)
  sys.exit(1)

if args.field_value is None:
  print('Error: Value to get the data is not specified', file=sys.stderr)
  sys.exit(1)

if args.input and args.input != '-':
    input_file = open(args.input)
else:
    input_file = sys.stdin

file_separator = args.input_field_sep
colname = args.colname
specified_value = args.field_value
outsep = args.output_field_sep

header_line = input_file.readline().rstrip("\n").split(file_separator)

if len(header_line) == 1:
  print("File has only one column, separator correct?", file=sys.stderr)

try:
  position = header_line.index(colname)
except ValueError:
  print("Error: Column name '%s' not found" % colname, file=sys.stderr)
  sys.exit(1)
else:
  if args.H == True:
    print(outsep.join(header_line))
  for line in input_file:
    line_list = line.rstrip("\n").split(file_separator)
    if line_list[position] == specified_value:
      outline = outsep.join(line_list)
      print(outline)
