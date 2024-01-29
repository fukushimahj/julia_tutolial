using Printf

# make initial condition
nmesh = 100   # mesh number
iend  = 200   # final step
ccfl  = 0.15
gamma = 5.0/3.0
iout  = 10     # file output

# label for physical value
M_rho   = 1
M_m     = 2
M_E     = 3

delx = 1.e0/nmesh
delt = 1.e-2
Ccfl = 0.1

folder = "./data/"
if !(isdir(folder))
  mkdir(folder)
end

function make_ini()
  
  v   = zeros(Float64, nmesh) # velocity 
  rho = zeros(Float64, nmesh) # 
  ein = zeros(Float64, nmesh) # internal energy

  rho0 = 1.e0
  p0   = 1.e0

  ## density
  for i = 1:Int(nmesh*0.5)
    rho[i] = rho0
  end
  for i = Int(nmesh*0.5):nmesh
    rho[i] = rho0*0.5
  end

  ## velocity
  for i=1:nmesh
    v[i] = 0.e0
  end

  for i=1:nmesh
    ein[i] = p0*(rho[i]/rho0)/(gamma-1.e0)
  end

  return v, rho, ein
end

function convert_v2u(v, rho, ein)
  
  U       = zeros(Float64, 3, nmesh)
  for i=1:nmesh
    U[M_rho,i] = rho[i]
    U[M_m,i]   = rho[i]*v[i]
    U[M_E,i]   = ein[i] + 0.5*rho[i]*v[i]^2
  end
  return U
end

function convert_u2v(U)

  v   = zeros(Float64, nmesh) # velocity 
  rho = zeros(Float64, nmesh) # 
  ein = zeros(Float64, nmesh) # internal energy

  for i=1:nmesh
    rho[i] = U[M_rho, i]
    v[i]   = U[M_m,i]/rho[i]
    ein[i] = U[M_E,i] - 0.5*rho[i]*v[i]^2
  end
  return v, rho, ein
  
end

function get_F(U)

  F      = zeros(Float64, 3, nmesh)
  for i=1:nmesh

    rho = U[M_rho,i]
    v   = U[M_m,  i]/rho
    E   = U[M_E,i]
    p   = (E-0.5*rho*v^2)*(gamma-1.e0)

    F[1,i] = U[M_m,i]
    F[2,i] = rho*v^2 + p
    F[3,i] = (E+p)*v
  end
  return F
end

function get_FHLL(U, F)

  FHLL = zeros(Float64, 3, nmesh)

  dt   = 1.e50

  for i=2:nmesh

    rho_L = U[M_rho,i-1]
    v_L   = U[M_m,  i-1]/rho_L
    E_L   = U[M_E,i-1]
    p_L   = (E_L-0.5*rho_L*v_L^2)*(gamma-1.e0)
    c_L   = sqrt(gamma*p_L/rho_L)
    
    rho_R = U[M_rho,i]
    v_R   = U[M_m,  i]/rho_R
    E_R   = U[M_E,i]
    p_R   = (E_R-0.5*rho_R*v_R^2)*(gamma-1.e0)
    c_R   = sqrt(gamma*p_R/rho_R)

    S_L   = min(v_L-c_L, v_R-c_R)
    S_R   = min(v_L+c_L, v_R+c_R)

    if(S_L > 0)
      for lab=1:3     
        FHLL[lab, i] = F[lab, i-1]
      end
    elseif(S_R >= 0)
      for lab=1:3
        FHLL[lab, i] = (S_R*F[lab, i-1] - S_L*F[lab, i] + S_L*S_R*(U[lab, i] - U[lab, i-1]))/(S_R-S_L)
      end
    else
      for lab=1:3     
        FHLL[lab, i] = F[lab, i]
      end
    end

    # get dt 
    dt = min( dt, Ccfl*delx/(c_R+abs(v_R)))
  end

  return FHLL, dt
end


function main()

  # make initial condition
  v, rho, ein = make_ini()
  
  # initialize 
  i_step  = 0
  time    = 0.e0
  F       = zeros(Float64, nmesh)
  F_HLL   = zeros(Float64, nmesh)

  Unew = zeros(Float64, 3, nmesh)

  # データ書き込み -----------------------------
  dt = 0.0
  filename = folder*"/data_0.dat"
  f = open(filename, "w")
  for i=2:nmesh-1
    @printf(f, "%e %e %e %e %e\n", delx*i, dt, rho[i], v[i], ein[i])
  end
  close(f)
  # ----------------------------------------------


  while(true)

    # timestep ----------
    i_step += 1
    if(i_step > iend)
      break
    end
    # -------------------

    # boundary condition ------------
    # free boundary 
    v[1]    = v[2]
    rho[1]  = rho[2]
    ein[1]  = ein[2]

    v[nmesh]   = v[nmesh-1]
    rho[nmesh] = rho[nmesh-1]
    ein[nmesh] = ein[nmesh-1]
    # -------------------------------

    U = convert_v2u(v, rho, ein)
    F = get_F(U)
    FHLL, dt = get_FHLL(U, F)

    for i=2:nmesh-1
      for lab=1:3
        U[lab, i] = U[lab, i] - dt/delx*(FHLL[lab, i+1]-FHLL[lab,i])
      end
    end

    v, rho, ein = convert_u2v(U)

    # データ書き込み -----------------------------
    if(i_step % iout == 0)
      filename = folder*"/data_"*string(i_step)*".dat"
      f = open(filename, "w")
      for i=2:nmesh-1
        @printf(f, "%e %e %e %e %e\n", delx*i, dt, rho[i], v[i], ein[i])
      end
      close(f)
    end
    # ----------------------------------------------
  end
end

main()
