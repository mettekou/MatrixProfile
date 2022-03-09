namespace MatrixProfile.Benchmarks

open System

open BenchmarkDotNet.Attributes
open BenchmarkDotNet.Jobs

open MatrixProfile

[<SimpleJob(RuntimeMoniker.Net50, baseline = true)>]
[<RPlotExporter>]
type VaryingWindowSize() =
    member val Jobs = 4
    
    [<Params(1024, 2048, 4096, 8192, 16384, 32768)>]
    member val M = 0 with get, set
    
    member val N = pown 2 15
    
    member val Series = Array.empty<float> with get, set
    
    [<GlobalSetup>]
    member this.Setup() =
        let random = Random()
        this.Series <- Array.init this.N (fun _ -> random.NextDouble())

    [<Benchmark>]
    member this.SelfJoin() =
        Mpx.selfJoin this.Series this.M true this.Jobs |> ignore