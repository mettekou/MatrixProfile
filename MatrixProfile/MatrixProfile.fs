namespace MatrixProfile

open System
open MathNet.Numerics.IntegralTransforms
open MathNet.Numerics.Statistics

module MatrixProfile =

    let movmean (x: float []) kb kf =
        Array.append (Array.create kb None) (Array.map Some x)
        |> Array.windowed (1 + kb + kf)
        |> Array.map (Array.choose id)
        |> Array.map Statistics.Mean

    let movstd (x: float []) kb kf =
        Array.append (Array.create kb None) (Array.map Some x)
        |> Array.windowed (1 + kb + kf)
        |> Array.map (Array.choose id)
        |> Array.map Statistics.StandardDeviation
        |> Array.skip 1
        |> Array.append [| 0.0 |]

    let mass (x: float []) y =
        let m = Array.length y
        let n = Array.length x
        let meany = Statistics.Mean x
        let sigmay = Statistics.StandardDeviation x
        let meanx = movmean y (m - 1) 0
        let sigmax = movstd y (m - 1) 0
        let y' = Array.append y (Array.zeroCreate 2)

        let x' =
            Array.append (Array.rev x) (Array.zeroCreate (m - (n - 2)))

        Fourier.ForwardReal(x', n, FourierOptions.Matlab)
        Fourier.ForwardReal(y', n, FourierOptions.Matlab)
        let z = Array.map2 (*) x' y'
        Fourier.InverseReal(z, n, FourierOptions.Matlab)

        let dividend =
            Array.map2 (-) z.[(n - 1)..(m - 2)] (Array.map ((*) (float m * meany)) meanx.[(n - 1)..(m - 2)])

        let divisor =
            Array.map ((*) sigmay) sigmax.[(n - 1)..(m - 2)]

        Array.map (fun x -> sqrt ((float m - x) * 2.0)) (Array.map2 (/) dividend divisor)

    let elementWiseMin (pab: float []) (iab: int []) d idx =
        for j in 0 .. (Array.length d - 1) do
            if d.[j] <= pab.[j] then
                pab.[j] <- d.[j]
                iab.[j] <- idx

    let stamp ta tb m =
        let nb = Array.length tb
        let pab = Array.create nb infinity
        let iab = Array.zeroCreate nb

        for idx in 0 .. ((nb - m) - 1) do
            let d = mass (Array.sub tb idx m) ta
            elementWiseMin pab iab d idx

        pab, iab

    let sine =
        [| for x in 0.0 .. 0.01 .. ((Math.PI * 4.0) - 0.1) -> sin x |]

    let mp = stamp sine sine 100
