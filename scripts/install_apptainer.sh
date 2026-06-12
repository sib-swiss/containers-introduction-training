sudo apt update
sudo apt install -y software-properties-common

sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer

apptainer_path=`which apptainer`
# Set up apptainer profile in app armor to prevent permission issues
sed -i "s|APPTAINER_PATH|$apptainer_path|" setup_scripts/apptainer_profile.txt # I use the pipe as separator because apptainer_path contains slashes, which messes up sed
sudo cp scripts/apptainer_profile.txt /etc/apparmor.d/apptainer
sudo systemctl reload apparmor
