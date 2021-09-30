# Rule that creates a log file with the files and other parameters specified in the config file
# Has date that the pipeline was ran

from datetime import date

running_date = date.today()
running_date = running_date.strftime("%d-%m-%Y")
files_log = running_date + "_files.txt"

rule create_file_log:
    output:
        files_log
    run:
        with open(output[0], "w") as outfile:
            outfile.write(f"Files used to run {pipeline} on {running_date}\n")
            for key, value in config.items():     
                outfile.write(f"{key}: {value}\n")
