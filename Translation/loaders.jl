using MAT

function loadInData()
    file = matopen("Data/simulationData.mat")

    data = read(file, "data")

    close(file)

    return data
end

function loadStfp()
    file = matopen("Data/stfp.mat")

    stfp = read(file,"stfp")

    close(file)

    return stfp
end