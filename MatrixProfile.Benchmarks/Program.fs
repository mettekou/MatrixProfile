namespace MatrixProfile.Benchmarks

open BenchmarkDotNet.Running

module Program =

    [<EntryPoint>]
    let main _ =
        BenchmarkRunner.Run<VaryingSeriesSize>() |> ignore
        BenchmarkRunner.Run<VaryingWindowSize>() |> ignore
        0