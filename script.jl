using Veronese

function vhs_number(k, d)
    println("$k, $d")
    vhs = d + 1
    while true
        print("   trying $vhs... ")
        if vhs_identifiable(k * d, k, d, vhs; verbose=true)
            println("  success")
            return vhs
        else
            println("  failed")
        end
        vhs += 1
    end 
end

function ladmc_number(k, d)
    ladmc = d + 1 
    while true
        if binomial(ladmc + 1, 2) > k * binomial(d + 1, 2)
            return ladmc
        end
        ladmc += 1
    end
end
