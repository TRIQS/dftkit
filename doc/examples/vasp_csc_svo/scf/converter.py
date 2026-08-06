from triqs_dftkit.vasp import Converter
import triqs_dftkit.vasp.plovasp.converter as plo_converter

# Generate and store PLOs
plo_converter.generate_and_output_as_text('plo.cfg', vasp_dir='./')

# run the converter
Converter = Converter(filename = 'vasp',proj_or_hk='proj')
Converter.convert_dft_input()
