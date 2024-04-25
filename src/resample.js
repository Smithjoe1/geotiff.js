/**
 * @module resample
 */

function copyNewSize(array, width, height, samplesPerPixel = 1) {
  return new (Object.getPrototypeOf(array).constructor)(width * height * samplesPerPixel);
}

/**
 * Resample the input arrays using nearest neighbor value selection.
 * @param {TypedArray[]} valueArrays The input arrays to resample
 * @param {number} inWidth The width of the input rasters
 * @param {number} inHeight The height of the input rasters
 * @param {number} outWidth The desired width of the output rasters
 * @param {number} outHeight The desired height of the output rasters
 * @returns {TypedArray[]} The resampled rasters
 */
export function resampleNearest(valueArrays, inWidth, inHeight, outWidth, outHeight) {
  const relX = inWidth / outWidth;
  const relY = inHeight / outHeight;
  return valueArrays.map((array) => {
    const newArray = copyNewSize(array, outWidth, outHeight);
    for (let y = 0; y < outHeight; ++y) {
      const cy = Math.min(Math.round(relY * y), inHeight - 1);
      for (let x = 0; x < outWidth; ++x) {
        const cx = Math.min(Math.round(relX * x), inWidth - 1);
        const value = array[(cy * inWidth) + cx];
        newArray[(y * outWidth) + x] = value;
      }
    }
    return newArray;
  });
}

// simple linear interpolation, code from:
// https://en.wikipedia.org/wiki/Linear_interpolation#Programming_language_support
function lerp(v0, v1, t) {
  return ((1 - t) * v0) + (t * v1);
}

/**
 * Resample the input arrays using bilinear interpolation.
 * @param {TypedArray[]} valueArrays The input arrays to resample
 * @param {number} inWidth The width of the input rasters
 * @param {number} inHeight The height of the input rasters
 * @param {number} outWidth The desired width of the output rasters
 * @param {number} outHeight The desired height of the output rasters
 * @returns {TypedArray[]} The resampled rasters
 */
export function resampleBilinear(valueArrays, inWidth, inHeight, outWidth, outHeight) {
  const relX = inWidth / outWidth;
  const relY = inHeight / outHeight;

  return valueArrays.map((array) => {
    const newArray = copyNewSize(array, outWidth, outHeight);
    for (let y = 0; y < outHeight; ++y) {
      const rawY = relY * y;

      const yl = Math.floor(rawY);
      const yh = Math.min(Math.ceil(rawY), (inHeight - 1));

      for (let x = 0; x < outWidth; ++x) {
        const rawX = relX * x;
        const tx = rawX % 1;

        const xl = Math.floor(rawX);
        const xh = Math.min(Math.ceil(rawX), (inWidth - 1));

        const ll = array[(yl * inWidth) + xl];
        const hl = array[(yl * inWidth) + xh];
        const lh = array[(yh * inWidth) + xl];
        const hh = array[(yh * inWidth) + xh];

        const value = lerp(
          lerp(ll, hl, tx),
          lerp(lh, hh, tx),
          rawY % 1,
        );
        newArray[(y * outWidth) + x] = value;
      }
    }
    return newArray;
  });
}


/**
 * Resample the input arrays using bicubic interpolation.
 * @param {TypedArray[]} valueArrays The input arrays to resample
 * @param {number} inWidth The width of the input rasters
 * @param {number} inHeight The height of the input rasters
 * @param {number} outWidth The desired width of the output rasters
 * @param {number} outHeight The desired height of the output rasters
 * @returns {TypedArray[]} The resampled rasters
 */
export function resampleBiCubic(valueArrays, inWidth, inHeight, outWidth, outHeight) {
  const relX = inWidth / outWidth;
  const relY = inHeight / outHeight;
  const newArrayLength = outWidth * outHeight;
  const newArray = new Array(valueArrays.length);

  // Function for bicubic interpolation
  function cubicInterpolate(p0, p1, p2, p3, x) {
    return 0.5 * (
      p0 * (-x + 2 * x * x - x * x * x) +
      p1 * (2 - 5 * x * x + 3 * x * x * x) +
      p2 * (x + 4 * x * x - 3 * x * x * x) +
      p3 * (-x * x + x * x * x)
    );
  }

  for (let i = 0; i < valueArrays.length; i++) {
    const array = valueArrays[i];
    const resultArray = new array.constructor(newArrayLength);

    for (let y = 0; y < outHeight; ++y) {
      for (let x = 0; x < outWidth; ++x) {
        const targetX = relX * x;
        const targetY = relY * y;

        const xFloor = Math.floor(targetX);
        const yFloor = Math.floor(targetY);

        const p = new Array(16);

        // Collect 16 neighboring pixels for bicubic interpolation
        for (let j = 0; j < 4; ++j) {
          for (let k = 0; k < 4; ++k) {
            const xIndex = xFloor - 1 + k;
            const yIndex = yFloor - 1 + j;

            let pixelValue;
            if (xIndex >= 0 && xIndex < inWidth && yIndex >= 0 && yIndex < inHeight) {
              pixelValue = array[yIndex * inWidth + xIndex];
            } else {
              pixelValue = 0; // If out of bounds, consider it as 0
            }

            p[j * 4 + k] = pixelValue;
          }
        }

        // Perform bicubic interpolation
        const xFraction = targetX - xFloor;
        const yFraction = targetY - yFloor;

        const interpolatedValue = cubicInterpolate(
          cubicInterpolate(p[0], p[1], p[2], p[3], xFraction),
          cubicInterpolate(p[4], p[5], p[6], p[7], xFraction),
          cubicInterpolate(p[8], p[9], p[10], p[11], xFraction),
          cubicInterpolate(p[12], p[13], p[14], p[15], xFraction),
          yFraction
        );

        resultArray[y * outWidth + x] = interpolatedValue;
      }
    }
    newArray[i] = resultArray;
  }

  return newArray;
}


/**
 * Resample the input arrays using the selected resampling method.
 * @param {TypedArray[]} valueArrays The input arrays to resample
 * @param {number} inWidth The width of the input rasters
 * @param {number} inHeight The height of the input rasters
 * @param {number} outWidth The desired width of the output rasters
 * @param {number} outHeight The desired height of the output rasters
 * @param {string} [method = 'nearest'] The desired resampling method
 * @returns {TypedArray[]} The resampled rasters
 */
export function resample(valueArrays, inWidth, inHeight, outWidth, outHeight, method = 'nearest') {
  switch (method.toLowerCase()) {
    case 'nearest':
      return resampleNearest(valueArrays, inWidth, inHeight, outWidth, outHeight);
    case 'bilinear':
    case 'linear':
      return resampleBilinear(valueArrays, inWidth, inHeight, outWidth, outHeight);
    case 'bicubic':
    case 'cubic':
      return resampleBiCubic( valueArrays, inWidth, inHeight, outWidth, outHeight);
    default:
      throw new Error(`Unsupported resampling method: '${method}'`);
  }
}

















/**
 * Resample the pixel interleaved input array using nearest neighbor value selection.
 * @param {TypedArray} valueArrays The input arrays to resample
 * @param {number} inWidth The width of the input rasters
 * @param {number} inHeight The height of the input rasters
 * @param {number} outWidth The desired width of the output rasters
 * @param {number} outHeight The desired height of the output rasters
 * @param {number} samples The number of samples per pixel for pixel
 *                         interleaved data
 * @returns {TypedArray} The resampled raster
 */
export function resampleNearestInterleaved(
  valueArray, inWidth, inHeight, outWidth, outHeight, samples) {
  const relX = inWidth / outWidth;
  const relY = inHeight / outHeight;

  const newArray = copyNewSize(valueArray, outWidth, outHeight, samples);
  for (let y = 0; y < outHeight; ++y) {
    const cy = Math.min(Math.round(relY * y), inHeight - 1);
    for (let x = 0; x < outWidth; ++x) {
      const cx = Math.min(Math.round(relX * x), inWidth - 1);
      for (let i = 0; i < samples; ++i) {
        const value = valueArray[(cy * inWidth * samples) + (cx * samples) + i];
        newArray[(y * outWidth * samples) + (x * samples) + i] = value;
      }
    }
  }
  return newArray;
}

/**
 * Resample the pixel interleaved input array using bilinear interpolation.
 * @param {TypedArray} valueArrays The input arrays to resample
 * @param {number} inWidth The width of the input rasters
 * @param {number} inHeight The height of the input rasters
 * @param {number} outWidth The desired width of the output rasters
 * @param {number} outHeight The desired height of the output rasters
 * @param {number} samples The number of samples per pixel for pixel
 *                         interleaved data
 * @returns {TypedArray} The resampled raster
 */
export function resampleBilinearInterleaved(
  valueArray, inWidth, inHeight, outWidth, outHeight, samples) {
  const relX = inWidth / outWidth;
  const relY = inHeight / outHeight;
  const newArray = copyNewSize(valueArray, outWidth, outHeight, samples);
  for (let y = 0; y < outHeight; ++y) {
    const rawY = relY * y;

    const yl = Math.floor(rawY);
    const yh = Math.min(Math.ceil(rawY), (inHeight - 1));

    for (let x = 0; x < outWidth; ++x) {
      const rawX = relX * x;
      const tx = rawX % 1;

      const xl = Math.floor(rawX);
      const xh = Math.min(Math.ceil(rawX), (inWidth - 1));

      for (let i = 0; i < samples; ++i) {
        const ll = valueArray[(yl * inWidth * samples) + (xl * samples) + i];
        const hl = valueArray[(yl * inWidth * samples) + (xh * samples) + i];
        const lh = valueArray[(yh * inWidth * samples) + (xl * samples) + i];
        const hh = valueArray[(yh * inWidth * samples) + (xh * samples) + i];

        const value = lerp(
          lerp(ll, hl, tx),
          lerp(lh, hh, tx),
          rawY % 1,
        );
        newArray[(y * outWidth * samples) + (x * samples) + i] = value;
      }
    }
  }
  return newArray;
}
/**
 * Resample the pixel interleaved input array using bicubic interpolation.
 * @param {TypedArray} valueArray The input array to resample
 * @param {number} inWidth The width of the input raster
 * @param {number} inHeight The height of the input raster
 * @param {number} outWidth The desired width of the output raster
 * @param {number} outHeight The desired height of the output raster
 * @param {number} samples The number of samples per pixel for pixel interleaved data
 * @returns {TypedArray} The resampled raster
 */
export function resampleBiCubicInterleaved(valueArray, inWidth, inHeight, outWidth, outHeight, samples) {
  const relX = inWidth / outWidth;
  const relY = inHeight / outHeight;
  const newArrayLength = outWidth * outHeight * samples;
  const newArray = new valueArray.constructor(newArrayLength);

  // Function for bicubic interpolation
  function cubicInterpolate(p0, p1, p2, p3, x) {
    return 0.5 * (
      p0 * (-x + 2 * x * x - x * x * x) +
      p1 * (2 - 5 * x * x + 3 * x * x * x) +
      p2 * (x + 4 * x * x - 3 * x * x * x) +
      p3 * (-x * x + x * x * x)
    );
  }

  for (let y = 0; y < outHeight; ++y) {
    for (let x = 0; x < outWidth; ++x) {
      const targetX = relX * x;
      const targetY = relY * y;

      const xFloor = Math.floor(targetX);
      const yFloor = Math.floor(targetY);

      const p = new Array(16);

      // Collect 16 neighboring pixels for bicubic interpolation
      for (let j = 0; j < 4; ++j) {
        for (let k = 0; k < 4; ++k) {
          const xIndex = xFloor - 1 + k;
          const yIndex = yFloor - 1 + j;

          for (let i = 0; i < samples; ++i) {
            let pixelValue;
            if (xIndex >= 0 && xIndex < inWidth && yIndex >= 0 && yIndex < inHeight) {
              pixelValue = valueArray[(yIndex * inWidth + xIndex) * samples + i];
            } else {
              pixelValue = 0; // If out of bounds, consider it as 0
            }

            p[j * 4 + k] = pixelValue;
          }
        }
      }

      // Perform bicubic interpolation
      const xFraction = targetX - xFloor;
      const yFraction = targetY - yFloor;

      for (let i = 0; i < samples; ++i) {
        const interpolatedValue = cubicInterpolate(
          cubicInterpolate(p[0 * samples + i], p[1 * samples + i], p[2 * samples + i], p[3 * samples + i], xFraction),
          cubicInterpolate(p[4 * samples + i], p[5 * samples + i], p[6 * samples + i], p[7 * samples + i], xFraction),
          cubicInterpolate(p[8 * samples + i], p[9 * samples + i], p[10 * samples + i], p[11 * samples + i], xFraction),
          cubicInterpolate(p[12 * samples + i], p[13 * samples + i], p[14 * samples + i], p[15 * samples + i], xFraction),
          yFraction
        );

        newArray[(y * outWidth * samples) + (x * samples) + i] = interpolatedValue;
      }
    }
  }
  return newArray;
}


/**
 * Resample the pixel interleaved input array using the selected resampling method.
 * @param {TypedArray} valueArray The input array to resample
 * @param {number} inWidth The width of the input rasters
 * @param {number} inHeight The height of the input rasters
 * @param {number} outWidth The desired width of the output rasters
 * @param {number} outHeight The desired height of the output rasters
 * @param {number} samples The number of samples per pixel for pixel
 *                                 interleaved data
 * @param {string} [method = 'nearest'] The desired resampling method
 * @returns {TypedArray} The resampled rasters
 */
export function resampleInterleaved(valueArray, inWidth, inHeight, outWidth, outHeight, samples, method = 'nearest') {
  switch (method.toLowerCase()) {
    case 'nearest':
      return resampleNearestInterleaved(
        valueArray, inWidth, inHeight, outWidth, outHeight, samples,
      );
    case 'bilinear':
    case 'linear':
      return resampleBilinearInterleaved(
        valueArray, inWidth, inHeight, outWidth, outHeight, samples,
      );
    case 'bicubic':
    case 'cubic':
      return resampleBiCubicInterleaved(
        valueArray, inWidth, inHeight, outWidth, outHeight, samples,
      );
    default:
      throw new Error(`Unsupported resampling method: '${method}'`);
  }
}
