struct MPIRStats
    rhist::Vector
    TH::DataType
    TL::DataType
    TFact::DataType
    Meth::String
end

function MPIRStats(TH = Float64, TL = Float32, TFact = Float32)
    if TFact == TH
        MPStats = MPIRStats([], TH, TL, TFact, "Heavy IR")
    else
        MPStats = MPIRStats([], TH, TL, TFact, "IR")
    end
end
