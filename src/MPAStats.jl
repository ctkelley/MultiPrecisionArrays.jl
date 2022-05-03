struct MPGStats
   rhist::Vector
   khist::Vector{Vector}
   ididhist::Vector{Bool}
   TH::DataType
   TL::DataType
   TFact::DataType
   MethIR::DataType
   Meth::String
end

struct MPIRStats
   rhist::Vector
   TH::DataType
   TL::DataType
   TFact::DataType
   MethIR::DataType
   Meth::String
end

function MPIRStats(TH=Float64, TL=Float32, TFact=Float32)
    if TFact == TH
    MPStats=MPIRStats([],TH, TL, TFact,MPHArray,"Hungry IR")
    else
    MPStats=MPIRStats([],TH, TL, TFact,MPArray,"IR")
    end
end


function MPGStats()
    TH=Float64
    TL=Float16
    TFact=Float64
    MPStats=MPGStats([],[],[],TH, TL, TFact, MPGArray,"IRGM")
end

