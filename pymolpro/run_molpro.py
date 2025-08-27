#!env python
from pymolpro.molpro_input import InputSpecification, schema
import jsonschema
import pymolpro
import pathlib
import os
import argparse
import tempfile
import time

if __name__ == '__main__':

    def follow(thefile):
        '''generator function that yields new lines in a file
        '''
        while True:
            line = thefile.readline()
            if not line:
                time.sleep(0.1)
                continue
            yield line


    parser = argparse.ArgumentParser(
        prog='Molpro'
    )
    parser.add_argument('input', nargs=1)
    parser.add_argument('project_name', nargs='?', default=None)
    args = parser.parse_args()
    spec = InputSpecification(specification=args.input[0])
    # print('spec',spec)
    try:
        jsonschema.validate(instance=dict(spec), schema=schema)
        input = spec.molpro_input()
    except Exception as ex:
        # print('not taking spec',ex)
        input = args.input[0]

    # print('input',input)
    # print('json of input',json.dumps(dict(InputSpecification(input))))
    # exit()

    with tempfile.TemporaryDirectory() as tmpdir:
        if args.project_name is not None:
            fully_qualified_project_directory = pathlib.Path(args.project_name).absolute()
        else:
            fully_qualified_project_directory = (pathlib.Path(tmpdir) / 'project')
        location = fully_qualified_project_directory.parent
        project_name = fully_qualified_project_directory.name
        project = pymolpro.Project(project_name, location=location)
        project.write_input(input)
        project.run(wait=False)
        while not os.path.exists(project.filename('out', '', 0)):
            time.sleep(0.1)
        with open(project.filename('out', '', run=0), 'r') as o:
            lines = follow(o)
            for line in lines:
                print(line.rstrip())
                if 'Molpro calculation terminated' in line or (not line and o.tell() >= os.fstat(o.fileno()).st_size):
                    break
