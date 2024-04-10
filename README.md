# Cappable-seq RNA-seq

## Overview

This repository contains files and analysis for RNA-seq derived from cappable seq. It also contains `ansible` playbooks necessary to bring up a linux instance with the required conda environment.

## Server initialization and automation

### Rationale

Server initialization is a time consuming process requiring many user inputs. First the server has to be updated to the latest (or desired) version of the OS, then required programs and dependencies have to be installed. This sometimes requires data transfer of environment files or credentials. Finally, any other required data has to be transferred over. If there are security credentials involved they have to be manually entered at each step. This is a long process that keeps the user at the keyboard and engaged. I'd rather automate it.

### Strategy- Ansible playbooks

`Ansible playbooks` are `.yml` configuration files used by the program `ansible` to automate the execution of commands and transfer of files on other machines remote from the localhost. They're designed to bring up a remote server and install any packages and configurations needed while using user-supplied security credentials.

### Required startup tasks

1. **Update the server.** The OS on the server is chosen when the instance is created. I usually choose the latest version of ubuntu due to its high compatibility. The included playbooks update and reboot the OS, essentially using `apt update` and `apt upgrade` followed by `reboot`.

2. **Install `Miniconda3`.** This is our package manager for python packages. The playbook will obtain the latest from the official URL and install it to the home directory.

3. **Create a conda environment from the provided file.** To increase reproducibility, you should specify all the required packages in an `environment.yml` file named after the environment. In this case, use the `clipseq.yml` file that `miniconda` will use to create an environment named `clipseq`. This file is version-blind but specifies the channel order to prevent incompatibilities.

### Security

This section is technical. More detail can be found in the online docs for `ansible`.

To control access to these servers, oracle uses an `ssh_key` instead of a password. Obtain this file during instance creation and store it securely on your local machine. Limit access to the *owner only* by running `chmod 400` after downloading the file. `Ansible` can locate the key by creating a `sensitive_data.yml` file that includes the key's path and the corresponding server IP. Encrypt this file using `ansible-vault` and call the encrypted file using `vars_files` in the `playbook`.

Store the password to decrypt the `sensitive_data.yml` file in a `.password_file.txt` on the local machine. Again, restrict access using `chmod 400`. Don't use git to back this up. Protect it as it is the key to decrypting server security credientials. When running the `playbook` from the command line, use the `--vault-password-file` argument to provide the password.

### Server inventory

To specify what servers our `playbooks` can be run on, there is a `server_inventory.yml` file that contains the names of the servers, their IP addresses, and their usernames. `Ansible` needs this information to log in along with the security credentials described above. That means that this file needs to be edited to provide the IP address to the specific instance that you want to connect to and the name that `ansible` can recognize it by in the `playbook`.

When running the `playbook` from the command line, point `ansible` to the `server_inventory.yml` by using the `-i` flag.

## Genome

Reads will align to both human and WSN genomes as we expect viral RNAs to be captured in cappable-seq. Therefore, I will create a custom combo of hg38 annotated from RefSeq and our lab's WSN genome with custom annotations.

### Procedure

This will be detailed in a jupyter notebook in the `/genome` folder. `STAR` and `bowtie2` are both able to create indices from this custom reference genome. Annotations were based on sequences obtained from GenBank. Annotation is done using the `gb2gtf.py` script and this procedure is detailed in the notebook.

## Experiments

### Organization

Experiments are organized into labeled subfolders in this directory but use the same reference genome in `/genome`. Each experiment contains its own reads and processing steps/qc reports
