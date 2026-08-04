#!/usr/bin/env python3
"""pick_cgc_genes.py

Simple helper to output the default gene list in config.yaml or a custom list.
"""
import sys
import yaml

if __name__ == '__main__':
    with open('config.yaml') as f:
        cfg = yaml.safe_load(f)
    for g in cfg.get('genes', []):
        print(g)
