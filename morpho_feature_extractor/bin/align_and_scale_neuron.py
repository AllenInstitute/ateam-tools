import sys
import swank
import argparse
from neuron_morphology import swc

def run(specimen_id, input_swc, output_swc, units=None):
    if units is None:
        units = swank.swc_units(input_swc)

    nrn = swc.read_swc(input_swc)
    aff = swank.align_and_scale_neuron(specimen_id, nrn, units)
    nrn.apply_affine(aff)
    nrn.save(output_swc)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--units', choices=['pixels','microns'], default=None)
    parser.add_argument('specimen_id', type=int)
    parser.add_argument('input_swc')
    parser.add_argument('output_swc')
    args = parser.parse_args()

    
    run(args.specimen_id, args.input_swc, args.output_swc, units=args.units)

if __name__ == "__main__": main()
