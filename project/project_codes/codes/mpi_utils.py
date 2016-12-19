BLOCK_LO = lambda proc_id,N,N_proc: proc_id*N // N_proc
BLOCK_HI = lambda proc_id,N,N_proc: BLOCK_LO(proc_id+1,N,N_proc)-1
BLOCK_SIZE = lambda proc_id,N,N_proc: \
        BLOCK_LO(proc_id+1,N,N_proc) - BLOCK_LO(proc_id,N,N_proc)
BLOCK_OWN = lambda k,N,N_proc: (N_proc*(k+1.0)-1) // N

