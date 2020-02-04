package shapeMixture

import breeze.linalg.{DenseMatrix, min}
import breeze.numerics.abs
import breeze.stats.distributions.RandBasis
import scala.util.control.Breaks._
import scala.annotation.tailrec
import scala.util.Sorting

object DTW {
  implicit val basis: RandBasis = RandBasis.withSeed(2)

  def getBestPathRec(costMatrix:DenseMatrix[Double]): List[(Int, Int)] = {

    var curRow = costMatrix.rows - 1
    var curCol = costMatrix.cols - 1
    var path = List((curRow, curCol))

    @tailrec
    def go(currentBestPath: List[(Int, Int)]): List[(Int, Int)] = {
      val curRow = currentBestPath.head._1
      val curCol = currentBestPath.head._2
      if(curCol==1 & curCol==1){
        currentBestPath
      } else {
        if(costMatrix(curRow-1,curCol-1) <= min(costMatrix(curRow,curCol-1),costMatrix(curRow-1,curCol) )){
          go((curRow-1, curCol-1) :: currentBestPath)

        } else if(costMatrix(curRow,curCol-1) <= costMatrix(curRow-1,curCol) ) {
          go((curRow, curCol-1) :: currentBestPath)
        } else {
          go((curRow-1, curCol) :: currentBestPath)
        }
      }
    }
    val res = go(path)
    res.map(e => (e._1-1,e._2-1))
  }

  def reduceByHalf(X: List[Double]):List[Double]={

    //    X With Duplicate Last ELement If Uneven Length
    val XCompleted:List[Double] = if(!(X.length%2==0)){
      X ++ List(X.last)
    } else {
      X
    }

    val indices = XCompleted.indices.filter(i => i%2==0)
    indices.map(i => (XCompleted(i)+XCompleted(i+1))/2D).toList

  }


  def dtw(X: List[Double], Y: List[Double]): (Double, List[(Int,Int)]) = {
    val n = X.length
    val m = Y.length

    var costMatrix = DenseMatrix.tabulate[Double](n+1,m+1){case (i,j) => Double.PositiveInfinity}
    costMatrix.update(0,0,0D)

    for(i <- 0 until n){
      for(j <- 0 until m){
        val dist = abs(X(i)- Y(j))
        val cost = dist + min(min(costMatrix(i,j+1),costMatrix(i+1,j)),costMatrix(i,j))
        costMatrix.update(i+1,j+1,cost)
      }
    }

    (costMatrix(n,m), getBestPathRec(costMatrix))
  }

  def DTWWithWindow(X: List[Double], Y: List[Double], window: List[(Int,Int)]= List[(Int,Int)]()): (Double, List[(Int,Int)]) = {
    val n = X.length
    val m = Y.length

    var costMatrix = DenseMatrix.tabulate[Double](n + 1, m + 1) { case (i, j) => Double.PositiveInfinity }
    costMatrix.update(0, 0, 0D)
    for ((i,j) <- window) {
      val dist = abs(X(i) - Y(j))
      val cost = dist + min(min(costMatrix(i, j + 1), costMatrix(i + 1, j)), costMatrix(i, j))
      costMatrix.update(i + 1, j + 1, cost)
    }

    (costMatrix(n, m), getBestPathRec(costMatrix))
  }

  def set(path: List[(Int,Int)])={
    val path2 = Sorting.stableSort(path, (e1: (Int, Int), e2: (Int, Int)) => e1._2 < e2._2).toList
    val path3 = Sorting.stableSort(path2, (e1: (Int, Int), e2: (Int, Int)) => e1._1 < e2._1)
    path3.toList.distinct
  }

  def expandWindow(path: List[(Int,Int)], len_x:Int, len_y:Int, radius:Int): List[(Int,Int)] = {
    var newPath = set(path)
    for ((i,j) <- path){
      val A = (-radius until radius +1).map(i+_)
      val B = (-radius until radius +1).map(j+_)
      for (a <- A) {
        for (b <- B) {
          newPath = newPath :+ (a,b)
        }
      }
    }

    newPath = set(newPath)
    var window_ = List[(Int,Int)]()

    for ((i,j) <- newPath) {
      window_ = window_ ++ List((2*i, 2*j)) :+
        (2*i ,2*j+1)   :+
        (2*i+1, 2*j)   :+
        (2*i+1, 2*j+1)
    }

    var window = List[(Int,Int)]()

    var start_j = 0
    for (i <- 0 until len_x){
      var new_start_j = -1
      breakable {
        for (j <- start_j until len_y){
          if(window_.contains((i, j))){
            window = window :+ (i, j)
            if (new_start_j == -1){
              new_start_j = j
            }
          } else if (new_start_j != -1) {
            break
          }}
        start_j = new_start_j
      }
    }

    window
  }

  def fastDTW(x: List[Double], y: List[Double], radius:Int): (Double, List[(Int,Int)]) = {

    val n = x.length
    val m = y.length

    val minTSsize = radius+2
    if(x.length<minTSsize | y.length<minTSsize) {
      DTW.dtw(x,y)
    } else {
      val shrunkX = reduceByHalf(x)
      val shrunkY = reduceByHalf(y)
      val lowResPath: List[(Int,Int)] = fastDTW(shrunkX,shrunkY,radius)._2
      val window = expandWindow(lowResPath,x.length,y.length,radius)
      DTWWithWindow(x,y,window)
    }
  }


}
