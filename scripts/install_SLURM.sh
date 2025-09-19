sudo apt update -y
sudo apt install slurmd slurmctld -y

sudo touch /etc/slurm/slurm.conf
# move the contents of slurm.conf to /etc/slurm/slurm.conf

sudo systemctl start slurmctld
sudo systemctl start slurmd

sudo scontrol update nodename=localhost state=idle

