namespace MatrixProfile.Benchmarks

open System

open BenchmarkDotNet.Attributes
open BenchmarkDotNet.Jobs

open MatrixProfile

[<SimpleJob(RuntimeMoniker.Net60, baseline = true)>]
[<RPlotExporter>]
type VaryingSeriesSize() =
    member val Jobs = 4
    
    member val M = pown 2 5
    
    [<Params(1024, 2048, 4096, 8192, 16_384, 32_768)>]
    member val N = 0 with get, set
    
    member val Series = Array.empty<float> with get, set
    
    [<GlobalSetup>]
    member this.Setup() =
        let random = Random()
        this.Series <- Array.init this.N (fun _ -> random.NextDouble())

    [<Benchmark>]
    member this.SelfJoin() =
        Mpx.selfJoin this.Series this.M true this.Jobs |> ignore


