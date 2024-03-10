struct MPIRStats
    rhist::Vector
    TH::DataType
    TF::DataType
    TFact::DataType
    Meth::String
end

function MPIRStats(TH = Float64, TF = Float32, TFact = Float32)
    if TFact == TH
        MPStats = MPIRStats([], TH, TF, TFact, "Heavy IR")
    else
        MPStats = MPIRStats([], TH, TF, TFact, "IR")
    end
end
