#! /bin/bash

config='local/config/'
OS='macOS' ## macOS 

echo $OS
if [ $OS == linux ]; then
    profile='/home/$USER/.bashrc'
elif [ $OS == macOS ]; then
    profile='/Users/$USER/.bash_profile'
fi
echo $profile
cat $config\direnv_allow.txt >> $profile
echo "export PATH=$PATH:./local/bin" > .envrc

cat venv_allow.txt >> $profile
echo "layout python-venv python3.5"  >> .envrc

cat $config\direnvrc_config.txt > ~/.config/direnv/direnvrc

source $profile
direnv allow ./

pip install -r $config\requirements.txt
