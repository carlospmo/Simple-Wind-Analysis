function PI = PowerInstantaneous(speed, Cp, radius, rho_air)

    PI = Cp * pi * radius^2 * rho_air * speed^3 / 2;
    
end 