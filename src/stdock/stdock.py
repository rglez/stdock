# Created by roy.gonzalez-aleman at 13/11/2023
import sys
import time

import stdock.config as cfg
import stdock.runner as runner


# todo: in mapdock-values, change the name of the section [docking-small] to [docking]
# todo: Put prints in the code to inform the user about the progress of the program
# todo: deprecate the support of filtered poses and scores in the program

def main():
    """ Main function of the program"""
    start = time.time()

    # =============================================================================
    # Processing the configuration file
    # =============================================================================
    if len(sys.argv) != 2:
        raise ValueError(
            '\nSTDOCK syntax is: python stdock.py path-to-config-file')
    config_path = sys.argv[1]
    # config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/map-from-values-then-dock_hur.cfg'

    params = cfg.allowed_parameters
    templates = cfg.allowed_templates
    config_args = cfg.STDConfig(config_path, params, templates).config_args
    print('Parameters have been correctly parsed')

    # =============================================================================
    # Running STDock
    # =============================================================================
    std_runner = runner.STDRunner(config_args)
    std_runner.run()
    print(f'Normal termination in {time.time() - start:.2f} secs')
