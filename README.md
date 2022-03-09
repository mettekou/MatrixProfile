# `MatrixProfile` [![](https://buildstats.info/nuget/MatrixProfile?includePreReleases=true)](https://www.nuget.org/packages/MatrixProfile)

`MatrixProfile` is a .NET Standard 2.0 library for time series analysis through the [matrix profile](https://www.cs.ucr.edu/~eamonn/MatrixProfile.html).
The name matrix profile refers to a data structure and the algorithms which produce it (currently only MPX in this library), developed by the Keogh and Mueen research groups at UC-Riverside and the University of New Mexico. The goal of this library is to make these algorithms accessible to both the novice and expert through standardization of core concepts, a simplistic API, and sensible default parameter values.

## Installation

### .NET CLI

To add `MatrixProfile` to a project using the .NET CLI, run the command:

```
dotnet add package MatrixProfile --version 0.1.0-alpha.23
```

Be sure to replace `0.1.0-alpha.23` by the version you want to install.

### Interactive

For interactive use, either in C# or F#, add the following directive to your script file or enter it into your REPL (e.g. `dotnet fsi`):

```fsharp
#r "nuget: MatrixProfile, 0.1.0-alpha.23"
```

Be sure to replace `0.1.0-alpha.23` by the version you want to use.

## Usage

The examples below are self-contained for interactive use, but are trivial to adjust for use in an application.

### From F#

```fsharp
#r "nuget: MatrixProfile, 0.1.0-alpha.23"

open System

open MatrixProfile

let random = Random()
let noisySeries = [| for x in 0. .. 0.01 .. 4. * Math.PI -> sin x + random.NextDouble() |]
let noisyQuery = [| for x in 0. .. 0.01 .. Math.PI -> sin x + random.NextDouble() |]

MatrixProfile.similarityJoin noisySeries noisyQuery (Array.length noisyQuery) true 2
```

### From C#

```csharp
#r "nuget: MatrixProfile, 0.1.0-alpha.23"

using System;

using MatrixProfile;

var random = new Random();
var noisySeries = Enumerable.Range(0, 1257).Select(i => Math.Sin(i * 0.01) + random.NextDouble()).ToArray();
var noisyQuery = Enumerable.Range(0, 315).Select(i => Math.Sin(i * 0.01) + random.NextDouble()).ToArray();

MatrixProfile.SimilarityJoin(noisySeries, noisyQuery, noisyQuery.Count(), true, 2);
```