# Test Fishing Problem
# Switched Nonlinear System
# with Uniform Switching Costs
# without Constraints
function fishing( ; solver=Bazinga.ZEROFPR( tol_optim=1e-4,
                                            max_iter=100,
                                            verbose=true,
                                            freq=1 ),
                        objtol = 1e-3,
                        primaltol = 1e-3 )

    println("Test: Fishing, with switching costs")

    # time interval
    t0 = 0.0
    tf = 12.0

    # control sequence
    uvec = Matrix( [repeat([0.0; 1.0], 4, 1); 0.0]' )

    # cost matrix
    C = [1.0 0.0 -1.0 0.0;
         0.0 1.0 0.0  -1.0]
    Q = Matrix( C' * C)

    # initial state
    x0 = [0.5; 0.7; 1; 1]

    # system dynamics
    function nldyn( x::Vector{Float64}, u::Vector{Float64} )
        n = length(x)
        f = zeros(n)
        f[1] = x[1] - x[1]*x[2] - 0.4*x[1]*u[1]
        f[2] = -x[2] + x[1]*x[2] - 0.2*x[2]*u[1]
        return f
    end

    function nldyn_deriv( x::Vector{Float64}, u::Vector{Float64} )
        df = [1.0-x[2]-0.4*u[1]       -x[1]                   0    0;
              x[2]                     -1+x[1]-0.2*u[1]       0    0;
              0                        0                      0    0;
              0                        0                      0    0]
        return df
    end

    # number of points on the fixed grid
    ngrid = 100

    # switching cost
    swc = 1e-3

    # STO problem with switching costs
    prob = scstoproblem(x0, nldyn, nldyn_deriv, uvec, ngrid=ngrid, t0=t0, tf=tf, Q=Q, swc=swc)

    out = solver( prob, L=1.0 )

    print( out )

    # Test Optimal Solution
    #=tauopt = gettau(m)
    @test isapprox(tauopt[1], 2.4434718158206916, atol=primaltol)
    @test isapprox(tauopt[2], 4.122362522221089, atol=primaltol)
    @test isapprox(tauopt[3], 4.433072844726073, atol=primaltol)
    @test isapprox(tauopt[4], 4.68252002334405,  atol=primaltol)
    @test isapprox(tauopt[5], 5.201667827750571, atol=primaltol)
    @test isapprox(tauopt[6], 5.369908894975762, atol=primaltol)
    @test isapprox(tauopt[7], 6.376845190917198, atol=primaltol)
    @test isapprox(tauopt[8], 6.47160340206027, atol=primaltol)

    # Test Optimal Value
    @test isapprox(getobjval(m), 1.3454355602054182, atol=objtol)=#

    println("Passed")
end
