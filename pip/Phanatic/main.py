def main():

    import os
    import sys
    import subprocess
    import argparse
    import random
    import pandas as pd

    # Did you know prompts
    prompts = [
        "Phanatic currently only supports short, paired end, illumina reads (default 'PE_illumina_150').",
        "I (J. Iszatt) made Phanatic to assemble S. aureus phage genomes as part of my PhD project.",
        "This is a short read assembler primarily designed for bacteriophage",
        "My favourite bacteriophage is a Silviavirus named Koomba-kaat_1",
        "You can use your own config file to customise the assembly and functions"
    ]
    random_prompt = random.choice(prompts)

    # Creating function to check directory path
    def valid_dir(dir_path):
        if not os.path.isdir(dir_path):
            raise argparse.ArgumentTypeError(
                f"{dir_path} is not a valid directory path")
        if not os.access(dir_path, os.R_OK):
            raise argparse.ArgumentTypeError(
                f"{dir_path} is not a readable directory")
        return dir_path

    def valid_file(file_path):
        if not os.path.isfile(file_path):
            raise argparse.ArgumentTypeError(
                f"{file_path} is not a valid file path")
        if not os.access(file_path, os.R_OK):
            raise argparse.ArgumentTypeError(
                f"{file_path} is not a readable file")
        return file_path

    # Parsing arguments
    image = 'iszatt/phanatic:2.2.2'
    parser = argparse.ArgumentParser(description=f"Easy short read assembly. Joshua J Iszatt: https://github.com/JoshuaIszatt")

    # Input/output options
    parser.add_argument('-i', '--input', type=valid_dir, help='Input reads files')
    parser.add_argument('-o', '--output', type=valid_dir, help='Direct output to this location')
    parser.add_argument('-r', '--reads', type=str, choices=['PE_illumina_150'], default='PE_illumina_150', help='Pipeline options')
    parser.add_argument('-c', '--config', type=valid_file, help='Use config file to customise assembly')
    parser.add_argument('-v', '--version', action="store_true", help='Print the docker image version')
    parser.add_argument('--show_console', action="store_true", help='Include this flag to write output to console')
    parser.add_argument('--manual', action="store_true", help='Enter container interactively')
    args = parser.parse_args()

    # Printing version
    if args.version:
        print(image)
        sys.exit(0)

    # Obtaining absolute paths if entered correctly
    if args.input and args.output:
        input_path = os.path.abspath(args.input)
        output_path = os.path.abspath(args.output)
    elif args.input and not args.output:
        sys.exit("No output directory specified\nRun --help to see options")
    elif args.output and not args.input:
        sys.exit("No input directory specified\nRun --help to see options")
    else:
        sys.exit("No input or output directories specified\nRun --help to see options")
        

    # Printing command variables
    print(
        f"Program run: {image}",
        f"Input path: {input_path}",
        f"Output path: {output_path}",
        f"Reads type: {args.reads}",
        ">>>",
        f"Did you know:",
        f"{random_prompt}",
        ">>>\n",
        sep='\n'
        )

    # Copying config file to output dir
    if args.config:
        print(f"Using {args.config} file \n")
        os.system(f"cp {args.config} {args.output}/config.ini")

    if args.show_console:
        docker = "docker run"
    else:
        docker = "docker run -d"

    # Running docker
    if args.manual: 
        os.system(f"docker exec -it \
            $(docker run -d \
            -v {input_path}:/assemble/input \
            -v {output_path}:/assemble/output \
            {image} sleep 1d) bash")
    else:
        command = ["%s -v %s:/assemble/input -v %s:/assemble/output %s /assemble/bin/assemble.sh" %
                (docker, input_path, output_path, image)]
        result = subprocess.Popen(command, shell=True)
        print(command)


if __name__ == "__main__":
    exit(main())