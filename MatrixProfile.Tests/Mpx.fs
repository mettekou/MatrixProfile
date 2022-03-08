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
        testCase "Small series similarity join multi-threaded" <| fun _ ->
            let ts = [| 1.; 2.; 3.; 1.; 2.; 3.; 4.; 5.; 6.; 0.; 0.; 1.; 1.; 2.; 2.; 4.; 5.; 1.; 1.; 9. |]
            let query = [| 0.; 0.; 1.; 1.; 2.; 2.; 4.; 5. |]
            let w = 4

            let desired = [|
                2.36387589e+00; 2.82842712e+00; 2.17957574e+00; 6.40728972e-01;
                6.40728972e-01; 6.40728972e-01; 3.26103392e+00; 3.61947699e+00;
                3.39984131e+00; 0.00000000e+00; 4.21468485e-08; 0.00000000e+00;
                4.21468485e-08; 0.00000000e+00; 2.82842712e+00; 3.57109342e+00;
                1.73771570e+00
            |]
            let desired_pi = [| 0; 1; 4; 1; 1; 1; 2; 1; 4; 2; 1; 2; 3; 4; 2; 1; 3 |]

            let profile = Mpx.similarityJoin ts query w true 4

            Array.iter2 (fun a e -> Expect.floatClose Accuracy.low a e "") profile.MatrixProfile desired
            Expect.sequenceEqual profile.MatrixProfileIndex desired_pi ""
      ]