if [ $(id -u) == 5658 ]; then
    export HOME=/home/default
else
    export HOME=/  
fi
export PATH=$PATH:/usr/lib64/mpich/bin 
