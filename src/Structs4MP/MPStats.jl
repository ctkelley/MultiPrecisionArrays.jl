struct MPIRStats
    rhist::Vector
    TW::DataType
    TF::DataType
    TFact::DataType
    Meth::String
end

function MPIRStats(TW = Float64, TF = Float32, TFact = Float32)
    if TFact == TW
        MPStats = MPIRStats([], TW, TF, TFact, "Heavy IR")
    else
        MPStats = MPIRStats([], TW, TF, TFact, "IR")
    end
end
