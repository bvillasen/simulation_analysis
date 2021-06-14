


def Get_Configuration( n_total ):
  nx, ny, nz = 1, 1, 1
  if n_total == 8:  nx, ny, nz = 2, 2, 2
  if n_total == 16: nx, ny, nz = 4, 2, 2
  if n_total == 32: nx, ny, nz = 4, 4, 2
  if n_total == 64: nx, ny, nz = 4, 4, 4
  if n_total == 128: nx, ny, nz = 8, 4, 4
  if n_total == 256: nx, ny, nz = 8, 8, 4
  if n_total == 512: nx, ny, nz = 8, 8, 8
  if n_total == 1024: nx, ny, nz = 16, 8, 8
  if n_total == 2048: nx, ny, nz = 16, 16, 8
  if n_total == 4096: nx, ny, nz = 16, 16, 16
  if n_total == 8192: nx, ny, nz = 32, 16, 16
  if n_total == 16384: nx, ny, nz = 32, 32, 16
  if nx * ny * nz != n_total: 
    print( f'ERROR: Domain mismatch  {n_total}: {nx} {ny} {nz} ')
    return None
  return ( nx, ny, nz )