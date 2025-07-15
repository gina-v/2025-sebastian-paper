#! /usr/bin/env python

import os
import sys
import defusedxml.ElementTree as DET
import csv
import argparse


def display_tree(element, level=0):
    """
    Display the XML tree for posterity

    Parameters:
    element (str): The element of the XML tree to display
    level (int): the level to begin display
    """
    indent = "  " * level
    print(f"{indent}<{element.tag} {element.attrib}>")
    for child in element:
        display_tree(child, level + 1)

def process_xml(xml_file, output_csv):
    """
    Process an XML file and append its extracted statistics to a CSV file.

    Parameters:
    xml_file (str): Path to the input XML file.
    output_csv (str): Path to the output CSV file.
    """
    with open(output_csv, mode='a', newline='') as file:
        writer = csv.writer(file, delimiter=',')

        try:
            tree = DET.parse(xml_file)
            root = tree.getroot()

            print(f"Root element: <{root.tag} {root.attrib}>")

            print("\nDisplaying XML tree")
            display_tree(root)

            run_accession = root.attrib.get('accession')
            spot_count = root.attrib.get('spot_count')
            base_count = root.attrib.get('base_count')

            if os.path.exists(output_csv):
                with open(output_csv, mode='r', newline='') as file:
                    reader = csv.reader(file)
                    existing_accessions = {row[1] for row in reader}
            else:
                existing_accessions = set()

            if run_accession not in existing_accessions:
                if run_accession and spot_count and base_count:

                    print(f"\nWriting {spot_count} reads and {base_count} bases to {run_accession}")
                    writer.writerow([xml_file, run_accession, spot_count, base_count])

                else:
                    sequence_table = root.find(".//Table[@name='SEQUENCE']/Statistics")
                    accession = root.attrib.get('accession')
                    read_count = sequence_table.find('Rows').attrib['count']
                    base_count = sequence_table.find('Elements').attrib['count']

                    print(f"\nWriting {read_count} reads and {base_count} bases to {accession}")
                    writer.writerow([xml_file, accession, read_count, base_count])
            else:
                print(f"\n{run_accession} already in {output_csv}. Skipping...")

        except Exception as e:
            print(f"Error processing {xml_file}: {e}")

def main():
    p = argparse.ArgumentParser()

    p.add_argument('xml', help='XML file to parse')
    p.add_argument('-o', '--output', help='Appended line of stats to output file')

    args = p.parse_args()

    process_xml(args.xml, args.output)

if __name__ == '__main__':
    sys.exit(main())
