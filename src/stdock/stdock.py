# Created by roy.gonzalez-aleman at 13/11/2023
import sys
import time

import config as cfg
import std

start = time.time()

# =============================================================================
# Processing the configuration file
# =============================================================================
# if len(sys.argv) != 2:
#     raise ValueError('\nSTDOCK syntax is: python stdock.py path-to-config-file')
# config_path = sys.argv[1]
config_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/config.cfg'
params = cfg.allowed_parameters
templates = cfg.allowed_templates
config_args = cfg.STDConfig(config_path, params, templates).config_args
print('Parameters have been correctly parsed')

# =============================================================================
# STD-NMR related processing
# =============================================================================
std_runner = std.STDRunner(config_args)

print(f'Normal termination in {time.time() - start:.2f} secs')

# %%
# =============================================================================
# Docking simulations
# =============================================================================
