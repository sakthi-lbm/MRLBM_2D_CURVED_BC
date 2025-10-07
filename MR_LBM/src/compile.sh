LT=D2Q9

# Set CompCap if it's not already defined
if [ -z "$CompCap" ]; then
    CompCap=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | tr -d '.')
    if [ -z "$CompCap" ]; then
        echo "Error: Unable to determine compute capability."
        exit 1
    fi
fi
rm -f ./../*sim_D2Q9_sm75
rm -r ./../CYLINDER/$1/
    nvcc -gencode arch=compute_${CompCap},code=sm_${CompCap} -rdc=true -O3 --restrict \
        *.cu \
        -lcudadevrt -lcurand -o ./../$1sim_${LT}_sm${CompCap}

cd ../
./$1sim_${LT}_sm${CompCap}
