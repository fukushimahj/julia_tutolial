using Printf

# make initial condition
nmesh = 100   # mesh number
iend  = 100   # final step
c     = 1.0   # speed of scaler
delt  = 1.e-1 # timestep
delx  = 1.e0  # cell width
iout  = 20     # file output

if !(isdir("./data/"))
  mkdir("./data/")
end

function make_ini()
  u = zeros(Float64, nmesh)

  for i = 1:Int(nmesh*0.5)
    u[i] = 1.e0
  end
  for i = Int(nmesh*0.5):nmesh
    u[i] = 0.5e0
  end
  return u
end 

function main()

  u = make_ini()
  i_step = 0
  time = 0.e0

  flux = zeros(Float64, nmesh)

  while(true)
    
    # timestep --------
    i_step += 1
    if(i_step > iend)
      break
    end
    # -----------------

    for i=1:nmesh-1
      flux[i] = 0.5*c*(u[i]+u[i+1]) # FTCS scheme
    end

    #boundary conditions
    u[1] = 1.e0
    u[nmesh] = 0.5

    for i=2:nmesh-1
      u[i] = u[i] - delt/delx*(flux[i]-flux[i-1])
    end
    time+= delt

    # データ書き込み -----------------------------
    if(i_step % iout == 0)
      filename = "./data/data_"*string(i_step)*".dat"
      f = open(filename, "w")
      for i=2:nmesh-1
        @printf(f, "%e %e\n", delx*i, u[i])
      end
      close(f)
    end
    # ----------------------------------------------



  end

end 
main()



