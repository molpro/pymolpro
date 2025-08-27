import json
import os
import re
import string

import jsonschema

from pymolpro import molpro_input
from pymolpro.molpro_input import equivalent, canonicalise, InputSpecification, _supported_methods
import pytest


@pytest.fixture
def methods():
    molpro_input._supported_methods = ['RHF', 'CCSD', 'RKS', 'CASSCF', 'MRCI', 'UHF', 'UKS', 'OCC', 'OPTG',
                                      'FREQUENCIES', 'THERMO']
    return molpro_input._supported_methods
    # yield supported_methods


def test_file(methods, tmpdir):
    test_file = tmpdir / 'test-molpro_input.inp'
    test_text = 'Geometry={F;H,F,1.7};geometry=hf.xyz;basis=cc-pVTZ !some comment;rhf\nccsd\n'
    with open(test_file, 'w') as f:
        f.write(test_text)
    assert InputSpecification(test_text).without_defaults == InputSpecification(test_file).without_defaults


def test_create_input(methods):
    for spec in [
        {'geometry': 'F\nH,F,1.7',
         'basis': {'default': 'cc-pVTZ', 'elements': {}},
         # 'steps': [{'command': 'rks', 'options': ['b3lyp']}, {'command': 'ccsd'}],
            'method': ['rks,b3lyp','ccsd'],
         'hamiltonian': 'AE',
         },
    ]:
        specification = InputSpecification(specification=spec)

        # print('initial specification',specification)
        # print('created input',specification.create_input(),'---')
        # print('new_specification', InputSpecification(specification.create_input()))
        assert InputSpecification(specification.molpro_input()).with_defaults == specification.with_defaults
        assert InputSpecification(specification.molpro_input()).without_defaults == specification.without_defaults

    for input in [
        'Geometry={F;H,F,1.7};basis={default=cc-pVTZ,h=cc-pVDZ} !some comment;{ks,b3lyp};{ccsd}\n',
        'Geometry={F;H,F,1.7};basis={default=cc-pVTZ,h=cc-pVDZ} !some comment;{ks,b3lyp};ccsd\n',
        'Geometry={\nF;H,F,1.7};basis={default=cc-pVTZ,h=cc-pVDZ} !some comment;{ks,b3lyp};ccsd\n',
        'Geometry={\nF;H,F,1.7\n};basis={default=cc-pVTZ,h=cc-pVDZ} !some comment;{ks,b3lyp};ccsd\n',
        'Geometry={\nF\nH,F,1.7\n};basis={default=cc-pVTZ,h=cc-pVDZ} !some comment;{ks,b3lyp};ccsd\n',
        'geometry={\nHe\n}\nhf',
        'geometry={\nHe\n}\nhf\nccsd',
        'geometry=thing.xyz;rhf',
        'geometry={H};uhf',
        'geometry={H};{uhf}',
        'geometry={H};{rhf}',
        'geometry={H};{hf}',
        'geometry={H};{uks,b3lyp};ccsd',
        'geometry={H};{rks,b3lyp};ccsd',
        'geometry={H};{ks,b3lyp};ccsd',
        'geometry={H};uks,b3lyp;ccsd',
        'geometry={H};uks,b3lyp',
        'geometry={H};rks,b3lyp',
        'geometry={H};ks,b3lyp',
        'geometry={H};ccsd,option1,option2=thing,,',
        'geometry={H};ks,b3lyp,option1,option2,,',
        # 'geometry={F;H,F,1.7};basis\ndefault=cc-pvtz\nBe=svp\nend\nrhf;ccsd'
    ]:
        # print('new one---\n',input)
        specification = InputSpecification(input)
        regenerated_input = specification.molpro_input()
        regenerated_specification = InputSpecification(regenerated_input)
        assert regenerated_specification.without_defaults == specification.without_defaults
        if not equivalent(regenerated_input, input):
            print('input', input)
            print('specification', specification)
            print('regenerated_specification', regenerated_specification)
            print('regenerated_input', regenerated_input)
            canonicalised_input = canonicalise(input)
            print('canonicalised input', canonicalised_input, type(canonicalised_input))
            canonicalised_regenerated_input = canonicalise(regenerated_input)
            print('canonicalised regenerated_input', canonicalised_regenerated_input,
                  type(canonicalised_regenerated_input))
            assert canonicalised_input == canonicalised_regenerated_input


def test_recreate_input(methods):
    for input in [
        'geometry={Ne};basis={sp,ne,cc-pvdz;c};{df-rhf}',
    ]:
        specification = InputSpecification(input)
        assert 'basis' not in specification
        assert 'variables' not in specification
    for input in [
        'geomtyp=xyz;geometry={Ne};basis={default=cc-pV(T+d)Z-PP};gprint,basis;{df-rhf}',
        'geometry={\nHe\n}\nhf',
        'geometry={\nHe\n}\nhf\nccsd',
        'geometry={\nHe\n}\nhf\nccsd\n\n',
        'geometry={He}\nhf\nccsd\n\n',
        '\ngeometry={\nB\nH B 2.2\n}\nocc,5,1,1,context=mcscf\nrhf\ncasscf\nmrci',
        'geometry={He};rks,b3lyp',
        'geometry={He};{rks,b3lyp}',
        'geometry=newnewnew.xyz\nbasis=cc-pVTZ-PP\nrhf',
        'geometry=wed.xyz\nbasis=cc-pVTZ-PP\nset,charge=1,spin=1,thing=whatsit\nxx=yy,p=q\nrhf',
        'geometry={Ne};proc doit;{rhf};ccsd;endproc;{frequencies,proc=doit;thermo,temp=298;another}',
        'symmetry,nosym;geometry=ChemSpider-937.xyz;basis=cc-pV(T+d)Z-PP;gexpec,qm;gthresh,energy=1e-14;{df-rhf}',
        'geometry={Ne};basis=cc-pV(T+d)Z-PP;gprint,basis;{df-rhf}',
        'geometry={Ne;He,Ne,2};basis={default=cc-pV(T+d)Z-PP,He=vdz(s)};gprint,basis;{df-rhf}',
        'geometry={Ne};basis=cc-pV(T+d)Z-PP;gexpec,sm;{df-rhf}',
    ]:
        specification = InputSpecification(input)
        regenerated_input = specification.molpro_input()
        regenerated_specification = InputSpecification(regenerated_input)
        if not equivalent(regenerated_input, input) or regenerated_specification != specification:
            print('specification', specification)
            print('regenerated_specification', regenerated_specification)
            print('input', input)
            print('regenerated_input', regenerated_input)
            print('canonicalised input', canonicalise(input))
            print('canonicalised regenerated_input', canonicalise(regenerated_input))
        assert canonicalise(regenerated_input) == canonicalise(input)
        assert regenerated_specification == specification


def test_variables(methods):
    test_text = 'spin=2,charge=1! comment\nset,occ=[3,1,1] ! comments\n;Geometry={F;H,F,1.7};basis={default=cc-pVTZ,h=cc-pVDZ}\n{ks,b3lyp}!some comment;ccsd\n'
    specification = InputSpecification(test_text)
    # print('original input', test_text)
    # print('parsed specification', specification)
    # print('recreated input', create_input(specification))
    # print('parsed recreated input', InputSpecification(create_input(specification)))
    assert InputSpecification(specification.molpro_input()) == specification
    assert specification['variables']['spin'] == '2'
    assert specification['variables']['occ'] == '[3,1,1]'


def test_too_complex(methods):
    for test_text in [
        # 'geometry=a.xyz;geometry=b.xyz',
        # 'geometry=b.xyz;hf;ccsd;hf',
        'geometry=c.xyz;hf;basis=cc-pvtz;ccsd',
    ]:
        assert InputSpecification(test_text) == {}


def test_canonicalise(methods):
    for given, expected in {
        'geometry={\nHe\n}': 'geometry={he}\n',
        'a\n\n\nb\n': '{a}\n{b}\n',
        'basis={\ndefault=cc-pVTZ,h=cc-pVDZ\n} !some comment': 'basis=cc-pvtz,h=cc-pvdz !some comment\n',
        'basis={\ndefault=cc-pVTZ,h=cc-pVDZ\n} !some comment': 'basis=cc-pvtz,h=cc-pvdz\n',
    }.items():
        assert canonicalise(given) == expected
    for test_text in [
        'geometry={He}\nhf',
    ]:
        assert equivalent(test_text, InputSpecification(test_text).molpro_input())


def test_basis_qualities(methods):
    for test, quality in {
        'basis=cc-pVDZ': 2,
        'basis={default=cc-pVDZ}': 2,
        'basis={default=cc-pVDZ,H=cc-pVDZ}': 2,
        'basis={default=cc-pVTZ,H=cc-pVDZ}': 0,
    }.items():
        assert InputSpecification(test).basis_quality == quality


def test_basis_variants(methods):
    for test, outcome in {
        'basis=cc-pVDZ': 'basis=cc-pVDZ',
        'basis,cc-pVDZ': 'basis=cc-pVDZ',
        'basis=default=cc-pVDZ': 'basis=cc-pVDZ',
        'basis={default=cc-pVDZ}': 'basis=cc-pVDZ',
        'basis={cc-pVDZ}': 'basis=cc-pVDZ',
        'basis,default=cc-pVDZ': 'basis=cc-pVDZ',
        'basis,cc-pVDZ,h=cc-pVDZ(s)': 'basis=cc-pVDZ,H=cc-pVDZ(s)',
        'basis,cc-pVDZ,zR=cc-pVDZ(s),h=cc-pVTZ': 'basis=cc-pVDZ,Zr=cc-pVDZ(s),H=cc-pVTZ',
        'basis={cc-pVDZ,zR=cc-pVDZ(s),h=cc-pVTZ}': 'basis=cc-pVDZ,Zr=cc-pVDZ(s),H=cc-pVTZ',
    }.items():
        assert outcome in InputSpecification(test).molpro_input()


def test_method(methods):
    for test, outcome in {
        # '': None,
        'rhf': 'rhf',
        'rhf;ccsd': 'ccsd',
        'rhf;ccsd;optg;frequencies': 'ccsd',
        'rhf;ccsd;{optg};frequencies': 'ccsd',
        'rhf;ccsd;mrci;{optg};frequencies': 'mrci',
    }.items():
        specification = InputSpecification(test)
        assert specification.method == outcome
        for method in ['rhf', 'ccsd', 'ks']:
            specification.method = method
            assert specification.method == method


def test_method_options(methods):
    for test, outcome in {
        '': 'hf',
        'rhf': 'rhf',
        'rhf;ccsd': 'ccsd',
        'rhf;ccsd;optg;frequencies': 'ccsd',
        'rhf;ccsd;{optg};frequencies': 'ccsd',
        'rhf;ccsd;mrci;{optg};frequencies': 'mrci',
    }.items():
        specification = InputSpecification(test)
        assert specification.method == outcome
        if specification.method is not None:
            options = {'option1': 'value1', 'option2': 'value2'}
            specification.method_options = options
            assert specification.method_options == [k+'='+v for k,v in options.items()]


def test_job_type(methods):
    for test, outcome in {
        '': 'SP',
        'rhf': 'SP',
        'rhf;ccsd': 'SP',
        'proc ansatz;rhf;ccsd;endproc;optg,proc=ansatz;frequencies,proc=ansatz': 'OPT+FREQ',
        'proc ansatz;rhf;ccsd;endproc;{optg,proc=ansatz};frequencies,proc=ansatz': 'OPT+FREQ',
        'proc ansatz;rhf;ccsd;mrci;endproc;{optg,proc=ansatz};frequencies,proc=ansatz': 'OPT+FREQ',
        'proc ansatz;rhf;ccsd;{optg};mrci;endproc;frequencies,proc=ansatz': 'OPT+FREQ', # TODO really?
    }.items():
        specification = InputSpecification(test)
        assert specification['job_type'] == outcome
        for method in ['rhf', 'ccsd', 'ks']:
            specification.method = method
            assert specification['job_type'] == outcome
        for jt in molpro_input._default_job_type_commands:
            specification['job_type'] = jt
            assert specification['job_type'] == jt


def test_density_functional(methods):
    for test, outcome in {
        '': None,
        'hf': None,
        'ks': 'LDA',
        'rks,b3lyp': 'B3LYP',
        'rks,b3lyp,a=b': 'B3LYP',
    }.items():
        specification = InputSpecification(test)
        assert specification.density_functional == outcome
        if specification.density_functional:
            specification.density_functional = 'PBE'
            assert specification.density_functional == 'PBE'


def test_open_shell_electrons(methods, tmpdir):
    for test, outcome in {
        'geometry={He}': 0,
        'geometry={Li}': 1,
        'geometry={1;;Li 0 0 0}': 1,
        'geometry={2;;Li 0 0 0;H 1 0 0}': 0,
        'geometry={2;;Be 0 0 0;H 1 0 0}': 1,
        'geometry={N}': 1,
        'geometry={C}': 0,
        'geometry={C;He,C,1}': 0,
    }.items():
        specification = InputSpecification(test)
        # print(specification)
        assert specification.open_shell_electrons == outcome
        open_shell_xyz_file = tmpdir / 'open_shell_electrons.xyz'
        with open(open_shell_xyz_file, 'w') as f:
            f.write(re.sub('.*{(.*)}.*', r'\1', test))
        specification = InputSpecification(re.sub('{.*}', str(open_shell_xyz_file), test))
        assert specification.open_shell_electrons == outcome
        os.remove(open_shell_xyz_file)


def test_json():
    from jsonschema import validate
    good_strings = [
        '{}',
        '{"geometry":"He", "method":"hf"}',
        '{"geometry":"He", "method":"hf", "basis": {"default":"cc-pvdz"}}',
        '{"geometry":"He", "method":"hf", "basis": {"default":"cc-pvdz", "elements":{}}}',
        '{"geometry":"He", "method":"hf", "basis": {"default":"cc-pvdz", "elements":{"Cu":"cc-pVTZ-PP","Zn":"cc-pVTZ"}}}',
        '{"geometry":"He", "method": "hf", "hamiltonian":"PP"}',
        '{"geometry":"He", "method": "hf", "hamiltonian":"AE"}',
        '{"geometry":"He", "method": "hf", "orientation":"mass"}',
        '{ "method": [ "rhf",  "mp2"], "symmetry": "none", "geometry": "PubChem-962.xyz", "basis": {"default": "cc-pV(T+d)Z", "elements": {}}, "density_fitting": true, "orbitals": ["ibo", "pipek", "nbo", "boys"], "hamiltonian": "AE"}',
    ]
    good_strings.append(json.dumps(dict(molpro_input.InputSpecification('geometry={He}'))))
    bad_strings = [
        '{"geometry":"He", "hamiltonian":"pp"}',
        '{"geometry":"He", "orientation":"Mass"}',
        '{"geometry":"He", "method":"hf", "basis": {"defaul":"cc-pvdz"}}',
        '{"geometry":"He", "method":"hf", "basis": {"default":"cc-pvdz", "elements":["Cu","cc-pVTZ-PP"]}}',
        '{"geometry":"He", "method":"hf", "basis": {"default":"cc-pvdz", "elements":{"Cu":1,"Zn":2}}}',
        '{ "method": [ "rhf",  "mp2"], "symmetry": "none", "geometry": "PubChem-962.xyz", "basis": {"default": "cc-pV(T+d)Z", "elements": {}}, "density_fitting": true, "orbitals": ["ibo", "pipek", "nbo", "ibo"], "hamiltonian": "AE"}',
    ]
    with open('molpro_input.json', 'r') as f:
        schema = json.load(f)
    for string in good_strings:
        obj = json.loads(string)
        # print('string', string)
        # print('obj',obj)
        validate(instance=obj, schema=schema)
        # print(molpro_input.InputSpecification(specification=obj).molpro_input())
    for string in bad_strings:
        with pytest.raises(jsonschema.exceptions.ValidationError) as excinfo:
            obj = json.loads(string)
            # print('string', string)
            # print('obj',obj)
            validate(instance=obj, schema=schema)
    # print('schema',schema['properties']['orientation']['enum'])

def test_keyval():
    tests = {
        'method=ccsd' : '{"method" :\n"ccsd"}',
        'geometry="F;H,F,1.732"':'{"geometry":"F;H,F,1.732"}',
        'geometry="F;H,F,1.732", method=ccsd"':'{"geometry":"F;H,F,1.732","method":"ccsd"}',
        'geometry="F;H,F,1.732", method=ccsd,basis=cc-pvdz':'{"geometry":"F;H,F,1.732", "method":"ccsd", "basis":{"default" : "cc-pvdz"}}',
        'geometry="F;H,F,1.732", method=ccsd,basis="cc-pvdz"':'{"geometry":"F;H,F,1.732", "method":"ccsd", "basis":{"default" : "cc-pvdz"}}',
        'geometry="F;H,F,1.732", method=ccsd,basis="cc-pvdz,Cu=cc-vtz"':'{"geometry":"F;H,F,1.732", "method":"ccsd", "basis":{"default" : "cc-pvdz","elements":{"Cu":"cc-vtz"}}}',
        'geometry="F;H,F,1.732", method=ccsd,basis="cc-pvdz,Ni=cc-pvqz,Cu=cc-vtz"':'{"geometry":"F;H,F,1.732", "method":"ccsd", "basis":{"default" : "cc-pvdz","elements":{"Ni":"cc-pvqz","Cu":"cc-vtz"}}}',
    }
    for test in tests:
        converted = molpro_input._convert_keyval_to_json(test)
        assert re.sub('\s+','',converted) == re.sub('\s+','',tests[test])
        jsonschema.validate(instance=json.loads(converted), schema=molpro_input.schema)
