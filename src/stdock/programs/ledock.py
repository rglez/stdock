# Created by roy.gonzalez-aleman at 05/01/2024
class LeDock(Program):
    def create_config(self):
        _, min_coords, max_coords = self.get_rec_axis(output_min_max=True)
        x_min, y_min, z_min = min_coords
        x_max, y_max, z_max = max_coords

        out_dir = join(self.out_dir, 'ledock')
        os.makedirs(out_dir, exist_ok=True)
        config_name = join(out_dir, 'config.conf')
        lig_list_name = join(out_dir, 'ligands.list')
        with open(config_name, 'wt') as conf:
            conf.write(f'Receptor\n')
            conf.write(f'{self.rec_path}\n')
            conf.write(f'RMSD\n')
            conf.write(f'{self.rmsd_tol}\n')
            conf.write(f'Binding pocket\n')
            conf.write(f'{x_min} {x_max}\n')
            conf.write(f'{y_min} {y_max}\n')
            conf.write(f'{z_min} {z_max}\n')
            conf.write(f'Number of binding poses\n')
            conf.write(f'{self.num_poses}\n')
            conf.write(f'Ligands list\n')
            conf.write(f'ligands.list\n')
        with open(lig_list_name, 'wt') as llist:
            llist.write(f'{self.lig_path}\n')
        return config_name

    def get_commands(self):
        return f'{self.exe} {self.create_config()}'
