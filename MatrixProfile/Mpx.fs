namespace MatrixProfile

open System.Threading
open System.Threading.Tasks

module Mpx =

    let muinvn a w =
        let n = Array.length a
        let mutable x, z, c, a1, a2, a3, mu_a = 0., 0., 0., 0., 0., 0., 0.
        let profile_len = n - w + 1
        let h = Array.zeroCreate n
        let r = Array.zeroCreate n
        let mu = Array.zeroCreate profile_len
        let sigma = Array.zeroCreate profile_len
        let mutable p = a[0]
        let mutable s = 0.
        for i in 1..(w - 1) do
            x <- p + a[i]
            z <- x - p
            s <- s + ((p - (x - z)) + (a[i] - z))
            p <- x
        
        mu[0] <- (p + s) / float w
        for i in w..(n - 1) do
            x <- p - a[i - w + 1]
            z <- x - p
            s <- s + ((p - (x - z)) - (a[i - w] + z))
            p <- x

            x <- p + a[i]
            z <- x - p
            s <- s + ((p - (x - z)) + (a[i] - z))
            p <- x

            mu[i - w + 1] <- (p + s) / float w
        
        for i in 0..(profile_len - 1) do
            for j in i..((i + w) - 1) do
                mu_a <- a[j] - mu[i]
                h[j] <- mu_a * mu_a

                c <- float (pown 2 27 + 1) * mu_a
                a1 <- (c - (c - mu_a))
                a2 <- (mu_a - a1)
                a3 <- a1 * a2
                r[j] <- a2 * a2 - (((h[j] - a1 * a1) - a3) - a3)

            p <- h[i]
            s <- r[i]
            for j in (i + 1)..((i + w) - 1) do
                x <- p + h[j]
                z <- x - p
                s <- s + (((p - (x - z)) + (h[j] - z)) + r[j])
                p <- x

            if p + s = 0 then
                sigma[i] <- 0.
            else
                sigma[i] <- 1. / sqrt(p + s)
        
        (mu, sigma)

    /// <summary>Computes the matrix profile without using the Fourier transform</summary>
    let selfJoin series windowSize euclideanDistance nJobs =
        let n = Array.length series
        // the original implementation allows the minlag to be manually set
        // here it is always w / 4 similar to SCRIMP++
        let minlag = int (ceil (float windowSize / 4.0))
        let profileLen = n - windowSize + 1

        let mu, sigma = muinvn series windowSize
        
        let df = Array.zeroCreate profileLen
        let dg = Array.zeroCreate profileLen
        let mp = Array.create profileLen -1.
        let mpi = Array.create profileLen -1
        
        let tmpMp = Array2D.create nJobs profileLen -1.
        let tmpMpi = Array2D.create nJobs profileLen -1
        
        // this is where we compute the diagonals and later the matrix profile
        df[0] <- 0.
        dg[0] <- 0.
        Parallel.ForEach (seq {windowSize..(n - 1)}, fun i ->
            df[i - windowSize + 1] <- (0.5 * (series[i] - series[i - windowSize]))
            dg[i - windowSize + 1] <- (series[i] - mu[i - windowSize + 1]) + (series[i - windowSize] - mu[i - windowSize])) |> ignore
        
        let mutable lastThread = 0
        use lastThreadSemaphore = new SemaphoreSlim(0, 1)
        lastThreadSemaphore.Release() |> ignore
        
        Parallel.ForEach (seq {minlag + 1..(profileLen - 1)},
                          ParallelOptions(MaxDegreeOfParallelism = nJobs),
                          (fun () ->
                              lastThreadSemaphore.Wait()
                              let lastThread' = lastThread
                              lastThread <- lastThread + 1
                              lastThreadSemaphore.Release() |> ignore
                              lastThread', 0., 0., 0),
                          (fun diag state index (threadnum, c, c_cmp, col) ->
                            let mutable c' = c
                            let mutable c_cmp' = c_cmp
                            let mutable col' = col
                            for i in diag..((diag + windowSize) - 1) do
                                c' <- c' + ((series[i] - mu[diag]) * (series[i-diag] - mu[0]))

                            for offset in 0..(n - windowSize - diag) do
                                col' <- offset + diag
                                c' <- c' + df[offset] * dg[col'] + df[col'] * dg[offset]
                                c_cmp' <- c' * sigma[offset] * sigma[col']
                                
                                // update the distance profile and profile index
                                if c_cmp' > tmpMp[threadnum, offset] then
                                    tmpMp[threadnum, offset] <- c_cmp'
                                    tmpMpi[threadnum, offset] <- col'
                                
                                if c_cmp' > tmpMp[threadnum, col'] then
                                    if c_cmp' > 1.0 then
                                        c_cmp' <- 1.0
                                    tmpMp[threadnum, col'] <- c_cmp'
                                    tmpMpi[threadnum, col'] <- offset
                            threadnum, c', c_cmp', col'),
                            fun _ -> ()) |> ignore
        
        // combine parallel results...
        for i in 0..(Array2D.length1 tmpMp - 1) do
            for j in 0..(Array2D.length2 tmpMp - 1) do
                if tmpMp[i,j] > mp[j] then
                    if tmpMp[i, j] > 1.0 then
                        mp[j] <- 1.0
                    else
                        mp[j] <- tmpMp[i, j]
                    mpi[j] <- tmpMpi[i, j]
        
        // convert normalized cross correlation to euclidean distance
        if euclideanDistance then
            for i in 0..(profileLen - 1) do
                mp[i] <- sqrt(2.0 * float windowSize * (1.0 - mp[i]))
        
        { MatrixProfile = mp; MatrixProfileIndex = mpi }

    (*let abJoin ts query w cross_correlation =
        let mutable i, j, k, mx, threadnum = 0, 0, 0, 0, 0
        let n = Array.length ts
        let qn = Array.length query
        let mutable cov_, corr_, eucdist, mxdist = 0., 0., 0., 0.
        let profile_len = n - w + 1
        let profile_lenb = qn - w + 1
        let (mua, siga) = muinvn ts w
        let (mub, sigb) = muinvn query w
        let diff_fa = Array.zeroCreate profile_len
        let diff_ga = Array.zeroCreate profile_len
        let diff_fb = Array.zeroCreate profile_lenb
        let diff_gb = Array.zeroCreate profile_lenb
        let mp = Array.create profile_len -1.
        let mpi = Array.create profile_len -1
        let mpb = Array.create profile_lenb -1.
        let mpib = Array.create profile_lenb -1

        diff_fa[0] <- 0.
        diff_ga[0] <- 0.
        Parallel.ForEach (seq {w..(n - 1)}, fun i ->
        diff_fa[i - w + 1] <- 0.5 * (ts[i] - ts[i - w])
        diff_ga[i - w + 1] <- (ts[i] - mua[i - w + 1]) + (ts[i - w] - mua[i - w])
        ) |> ignore

        diff_fb[0] <- 0.
        diff_gb[0] <- 0.
        Parallel.ForEach (seq {w..(qn - 1)}, fun i ->
        diff_fb[i - w + 1] <- 0.5 * (query[i] - query[i - w])
        diff_gb[i - w + 1] <- (query[i] - mub[i - w + 1]) + (query[i - w] - mub[i - w])
        ) |> ignore

        Parallel.ForEach ((seq {0..profile_len}), (fun () -> { Tmp_Mp = Array.create profile_len -1.; Tmp_Mpi = Array.create profile_len -1; Tmp_Mpb = Array.create profile_lenb -1.; Tmp_Mpbi = Array.create profile_lenb -1}), (fun i state index { Tmp_Mp = tmp_mp; Tmp_Mpi = tmp_mpi; Tmp_Mpb = tmp_mpb; Tmp_Mpbi = tmp_mpib} ->
            let mx = if profile_len - i < profile_lenb then profile_len - i else profile_lenb
            cov_ <- 0.
            for j in i..((i + w) - 1) do
                cov_ <- cov_ + ((ts[j] - mua[i]) * (query[j-i] - mub[0]))
        for i in 0..(Array2D.length1 tmp_mp - 1) do
            for j in 0..(Array2D.length2 tmp_mp) do
                if tmp_mp[i,j] > mp[j] then
                    if tmp_mp[i, j] > 1.0 then
                        mp[j] <- 1.0
                    else
                        mp[j] <- tmp_mp[i, j]
                    mpi[j] <- tmp_mpi[i, j]
            for j in 0 .. (mx - 1) do
                k <- j + i
                cov_ <- cov_ + diff_fa[k] * diff_gb[j] + diff_ga[k] * diff_fb[j]
                corr_ <- cov_ * siga[k] * sigb[j]

                if corr_ > tmp_mp[k] then
                    tmp_mp[k] <- corr_
                    tmp_mpi[k] <- j

                if corr_ > tmp_mpb[j] then
                    tmp_mpb[j] <- corr_
                    tmp_mpib[j] <- k
            { Tmp_Mp = tmp_mp; Tmp_Mpi = tmp_mpi; Tmp_Mpb = tmp_mpb; Tmp_Mpbi = tmp_mpib}), fun _ -> ()) |> ignore

        Parallel.ForEach ((seq {0..profile_lenb}), (fun () -> { Tmp_Mp = Array.create profile_len -1.; Tmp_Mpi = Array.create profile_len -1; Tmp_Mpb = Array.create profile_lenb -1.; Tmp_Mpbi = Array.create profile_lenb -1}), (fun i state index { Tmp_Mp = tmp_mp; Tmp_Mpi = tmp_mpi; Tmp_Mpb = tmp_mpb; Tmp_Mpbi = tmp_mpib} ->
            let mx = if profile_lenb - i < profile_len then profile_lenb - i else profile_len
            cov_ <- 0.
            for j in i..((i + w) - 1) do
                cov_ <- cov_ + ((query[j] - mub[i]) * (ts[j-i] - mua[0]))

            for j in 0 .. (mx - 1) do
                k <- j + i
                cov_ <- cov_ + diff_fb[k] * diff_ga[j] + diff_gb[k] * diff_fa[j]
                corr_ <- cov_ * sigb[k] * siga[j]

                if corr_ > tmp_mpb[k] then
                    tmp_mpb[k] <- corr_
                    tmp_mpib[k] <- j

                if corr_ > tmp_mp[j] then
                    tmp_mp[j] <- corr_
                    tmp_mpi[j] <- k
            { Tmp_Mp = tmp_mp; Tmp_Mpi = tmp_mpi; Tmp_Mpb = tmp_mpb; Tmp_Mpbi = tmp_mpib}), fun _ -> ()) |> ignore

        for i in 0..(Array.length tmp_mp - 1) do
                for j in 0..range(tmp_mp.shape[1]) do
                    if tmp_mp[i,j] > mp[j] then
                        if tmp_mp[i, j] > 1.0 then
                            mp[j] = 1.0
                        else
                            mp[j] = tmp_mp[i, j]
                        
                        mpi[j] = tmp_mpi[i, j]
            
        for i in 0..(Array.length tmp_mpb - 1) do
                for j in 0..range(tmp_mpb.shape[1]) do
                    if tmp_mpb[i,j] > mpb[j] then
                        if tmp_mpb[i, j] > 1.0 then
                            mpb[j] = 1.0
                        else
                            mpb[j] = tmp_mpb[i, j]
                        
                        mpib[j] = tmp_mpib[i, j]
        (mp, mpi, mpb, mpib)*)