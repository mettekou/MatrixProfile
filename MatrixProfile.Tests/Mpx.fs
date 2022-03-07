namespace MatrixProfile.Tests

open Expecto

open MatrixProfile

module Mpx =

    [<Tests>]
    let tests =
      testList "MPX" [
        testCase "Small series self join Euclidean multi-threaded" <| fun _ ->
            let series = [| 0.; 1.; 1.; 1.; 0.; 0.; 2.; 1.; 0.; 0.; 2.; 1. |]
            let windowSize = 4
            let expected = [| 1.9550; 1.9550; 0.8739; 0.; 0.; 1.9550; 0.8739; 0.; 0. |]
            let expectedi = [| 4; 5; 6; 7; 8; 1; 2; 3; 4 |]
            let { MatrixProfile = actual; MatrixProfileIndex = actuali } = Mpx.selfJoin series windowSize true 8
            Array.iter2 (fun a e -> Expect.floatClose Accuracy.low a e "") actual expected
            Expect.sequenceEqual actuali expectedi ""
        testCase "Small series self join Pearson multi-threaded" <| fun _ ->
            let series = [| 0.; 1.; 1.; 1.; 0.; 0.; 2.; 1.; 0.; 0.; 2.; 1. |]
            let windowSize = 4
            let expected = [| 0.522232967867094; 0.522232967867094; 0.904534033733291; 1.; 1.; 0.522232967867094; 0.904534033733291; 1.; 1. |]
            let expectedi = [| 4; 5; 6; 7; 8; 1; 2; 3; 4 |]
            let { MatrixProfile = actual; MatrixProfileIndex = mpi } = Mpx.selfJoin series windowSize false 8
            Array.iter2 (fun a e -> Expect.floatClose Accuracy.low a e "") actual expected
            Expect.sequenceEqual mpi expectedi ""
      ]