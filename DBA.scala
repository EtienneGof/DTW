package shapeMixture

import breeze.linalg.{DenseMatrix, DenseVector, min, sum}
import breeze.numerics.{abs, pow}
import breeze.stats.distributions.RandBasis
import shapeMixture.DTW._
import scala.annotation.tailrec
import scala.collection.immutable.ListMap
import scala.util.Random

object DBA {
  implicit val basis: RandBasis = RandBasis.withSeed(2)

  def sampleDV(probabilities: DenseVector[Double]): Int = {
    val dist = (0 until probabilities.length) zip probabilities.toArray
    val threshold = scala.util.Random.nextDouble
    val iterator = dist.iterator
    var accumulator = 0.0
    while (iterator.hasNext) {
      val (cluster, clusterProb) = iterator.next
      accumulator += clusterProb
      if (accumulator >= threshold)
        return cluster
    }
    sys.error("Error")
  }

  def inertia(sequences: List[List[Double]], seq: List[Double]): Double = {
    sequences.map(element => DTW.dtw(element,seq)._1).map(pow(_,2)).sum
  }

  def medoid(data: List[List[Double]], sampleRatio: Double=1): List[Double] = {
    require(sampleRatio > 0 & sampleRatio <=1,"sampleRatio argument should be in )0,1]")
    println("Initialization: computing the inertia of each ts..")
    val dataSampled: List[List[Double]] = if(sampleRatio!=1){
      Random.shuffle(data).take((data.length * sampleRatio).toInt)
    } else {
      data
    }
    val inertiaPerTS = dataSampled.map(ts => inertia(data, ts))
    val idxMedoid: Int = inertiaPerTS.indices.minBy(inertiaPerTS)
    dataSampled(idxMedoid)
  }

  def updateDBA(sequences: List[List[Double]], center:List[Double]) : List[Double] = {
    val res: List[List[(Int,List[Double])]] = sequences.par.map(ts => {
      val resDTW = DTW.dtw(ts,center)
      resDTW._2.par.map(indices => (ts(indices._1),indices._2)).toList.groupBy(_._2).mapValues(xs => xs.map(_._1)).toList
    }).toList

    val resReduce  = res.reduce(_++_)
    val resReduceRegularized = resReduce.map(L => (L._1,List(L._2.sum/L._2.length.toDouble)))
    val resReduceBygroup = ListMap(resReduceRegularized.groupBy(_._1).toSeq.sortWith(_._1 < _._1):_*).mapValues(xs => xs.map(_._2).reduce(_++_)).values.toList
    resReduceBygroup.map(valList => sum(valList)/valList.length.toDouble)
  }

  def launch(sequences: List[List[Double]],
             initCenter: List[Double],
             nIter: Int =10,
             verbose:Boolean =false): List[Double] = {

    @tailrec
    def go(center: List[Double], currentIter: Int): List[Double] = {
      val newCenter = updateDBA(sequences,center)
      if(currentIter==nIter || center.equals(newCenter) ) {
        newCenter
      } else {
        go(newCenter, currentIter+1)
      }
    }
    go(initCenter, 1)
  }

}
