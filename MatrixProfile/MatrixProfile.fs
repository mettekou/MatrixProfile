namespace MatrixProfile

module MatrixProfile =
    
    let selfJoin = Mpx.selfJoin
    
    let similarityJoin = Mpx.similarityJoin
    
type MatrixProfile =
    static member SelfJoin(series, windowSize, euclideanDistance, nJobs) = MatrixProfile.selfJoin series windowSize euclideanDistance nJobs
    
    static member SimilarityJoin(series, query, windowSize, euclideanDistance, nJobs) = MatrixProfile.similarityJoin series query windowSize euclideanDistance nJobs