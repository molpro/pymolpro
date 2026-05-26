import numpy as np


class CubeData:
    def __init__(self, source: str | dict = None):
        if source is None:
            pass
        elif type(source) == str and source.endswith('.cube'):
            self.load_from_cube_file(source)
        elif isinstance(source, dict):
            for attr in ['title', 'natoms', 'origin', 'orbitals', 'nval', 'dimensions', 'cells', 'atoms', 'data',
                         'orbital_identifiers']:
                if attr in source:
                    setattr(self, attr, source[attr])
        else:
            raise ValueError('Invalid cube data source')

    def load_from_cube_file(self, filename: str):
        with open(filename, 'r') as f:
            self.title = [f.readline().strip() for i in range(2)]
            line = f.readline().strip().split()
            self.natoms = int(line[0])
            self.origin = tuple([float(line[k + 1]) for k in range(3)])
            self.orbitals = self.natoms < 0
            self.natoms = abs(self.natoms)
            if len(line) > 4:
                self.nval = int(line[4])
                if self.nval > 1:
                    raise ValueError('Multiple values not supported')
            self.dimensions = []
            self.cells = np.ndarray(shape=(3, 3), dtype=np.float64)
            for i in range(3):
                line = f.readline().strip().split()
                self.dimensions.append(int(line[0]))
                for j in range(3):
                    self.cells[i, j] = float(line[j + 1])
            self.atoms = []
            for atom in range(self.natoms):
                line = f.readline().strip().split()
                self.atoms.append({'atomic_number': int(line[0]), 'charge': float(line[1]),
                                   'xyz': tuple([float(line[k + 2]) for k in range(3)])})
            if self.orbitals:
                line = f.readline().strip().split()
                self.orbital_identifiers = [int(line[k + 1]) for k in range(int(line[0]))]
                if len(self.orbital_identifiers) > 1:
                    raise ValueError('Multiple orbitals not supported')
            self.data = np.ndarray(shape=(self.dimensions), dtype=np.float64)
            for i in range(self.dimensions[0]):
                for j in range(self.dimensions[1]):
                    k0 = 0
                    while True:
                        line = f.readline().strip().split()
                        for k, value in enumerate(line):
                            self.data[i, j, k + k0] = float(line[k])
                        k0 += len(line)
                        if k0 >= self.dimensions[2]:
                            break

    def __str__(self):
        buffer = self.title[0] + '\n' + self.title[1] + '\n'
        buffer += f'{self.natoms * (-1 if self.orbitals else 1)} {self.origin[0]} {self.origin[1]} {self.origin[2]} 1\n'
        for i in range(3):
            buffer += f'{self.dimensions[i]}'
            for j in range(3):
                buffer += f' {self.cells[i, j]}'
            buffer += '\n'
        for atom in self.atoms:
            buffer += f'{atom["atomic_number"]} {atom["charge"]}'
            for coord in atom["xyz"]:
                buffer += f' {coord}'
            buffer += '\n'
        if self.orbitals:
            buffer += f'{len(self.orbital_identifiers)}'
            for identifier in self.orbital_identifiers:
                try:
                    buffer += f' {int(float(identifier) * 10)}'
                except ValueError:
                    buffer += f' 1'
            buffer += '\n'
            for i in range(self.dimensions[0]):
                for j in range(self.dimensions[1]):
                    k0 = 0
                    while k0 < self.dimensions[2]:
                        for k in range(min(6, self.dimensions[2] - k0)):
                            buffer += f' {self.data[i, j, k0]:.5e}'.replace('e', 'E')
                            k0 += 1
                        buffer += '\n'

        return buffer

    def __len__(self):
        return self.dimensions[0] * self.dimensions[1] * self.dimensions[2]

    def bounds(self, threshold: float = 1e-4):
        if not isinstance(threshold, float):
            return None
        for i0 in range(self.dimensions[0]):
            if np.any(np.abs(self.data[i0, :, :]) >= threshold):
                break
        for i1 in reversed(range(self.dimensions[0])):
            if np.any(np.abs(self.data[i1, :, :]) >= threshold):
                break
        for j0 in range(self.dimensions[1]):
            if np.any(np.abs(self.data[:, j0, :]) >= threshold):
                break
        for j1 in reversed(range(self.dimensions[1])):
            if np.any(np.abs(self.data[:, j1, :]) >= threshold):
                break
        for k0 in range(self.dimensions[2]):
            if np.any(np.abs(self.data[:, :, k0]) >= threshold):
                break
        for k1 in reversed(range(self.dimensions[2])):
            if np.any(self.data[:, :, k1] >= threshold):
                break
        # print(self.dimensions)
        # print(i0, i1, j0, j1, k0, k1)
        return (
            (self.origin[0] + (i0-1) * self.cells[0, 0], self.origin[0] + (i1+1) * self.cells[0, 0]),
            (self.origin[1] + (j0-1) * self.cells[1, 1], self.origin[1] + (j1+1) * self.cells[1, 1]),
            (self.origin[2] + (k0-1) * self.cells[2, 2], self.origin[2] + (k1+1) * self.cells[2, 2])
        )
