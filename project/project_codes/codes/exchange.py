from mpi4py import MPI

def assign_aliases(u):
    return (u[1,1:-1], u[-2,1:-1], u[0,1:-1], u[-1,1:-1])

def exchange1d(cart1d, u, left, right):
    u_left, u_right, u_left_ghost, u_right_ghost = assign_aliases(u)
    cart1d.Send( u_left, dest=left )
    cart1d.Recv( u_right_ghost, source=right )
    cart1d.Send( u_right, dest=right )
    cart1d.Recv( u_left_ghost, source=left )
    return None

def exchange1d_paired(cart1d, u, left, right):
    pid = cart1d.Get_rank()
    u_left, u_right, u_left_ghost, u_right_ghost = assign_aliases(u)
    if pid%2 == 0:
        cart1d.Send( u_left, dest=left )
        cart1d.Recv( u_right_ghost, source=right )
        cart1d.Send( u_right, dest=right )
        cart1d.Recv( u_left_ghost, source=left )
    else:
        cart1d.Recv( u_right_ghost, source=right )
        cart1d.Send( u_left, dest=left )
        cart1d.Recv( u_left_ghost, source=left )
        cart1d.Send( u_right, dest=right )
    return None

def exchange1d_combined(cart1d, u, left, right):
    u_left, u_right, u_left_ghost, u_right_ghost = assign_aliases(u)
    cart1d.Sendrecv( sendbuf=u_left, dest=left, \
                     recvbuf=u_right_ghost, source=right )
    cart1d.Sendrecv( sendbuf=u_right, dest=right, \
                     recvbuf=u_left_ghost, source=left )
    return None

def exchange1d_nonblocking(cart1d, u, left, right):
    req = [ MPI.Request() for k in range(4) ]
    u_left, u_right, u_left_ghost, u_right_ghost = assign_aliases(u)
    req[0] = cart1d.Irecv( u_right_ghost, source=right )
    req[1] = cart1d.Irecv( u_left_ghost, source=left )
    req[2] = cart1d.Isend( u_left, dest=left )
    req[3] = cart1d.Isend( u_right, dest=right )
    return req
