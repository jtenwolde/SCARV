from scars import scars_queries, scars_codata, scars_assess, scars_diagnostics


def create_genomics_cnn():
    import keras
    import keras_genomics

    model = keras.models.Sequential()

    model.add(keras_genomics.layers.RevCompConv1D(filters=250, kernel_size=5,
                input_shape=(13,4), padding="same", activation="relu", kernel_initializer="he_normal"))

    model.add(keras_genomics.layers.core.DenseAfterRevcompConv1D(units=100, activation="relu"))

    return model



def create_combo(codata_shape):
    import keras
    from keras.models import Model
    from keras.layers import Input, Dense, concatenate

    seq_cnn = create_genomics_cnn()
    codata_input = Input(shape=codata_shape)

    merged = keras.layers.concatenate([seq_cnn.output, codata_input])

    next_layer = Dense(units=100, activation="relu")(merged)
    final_model_output = Dense(units=4, activation="softmax")(next_layer)

    final_model = Model(inputs=[seq_cnn.input, codata_input], outputs=final_model_output)

    return final_model



def create_genomics_cnn_seqonly():
    import keras
    import keras_genomics

    model = keras.models.Sequential()

    model.add(keras_genomics.layers.RevCompConv1D(filters=250, kernel_size=5,
                input_shape=(13,4), padding="same", activation="relu", kernel_initializer="he_normal"))

    model.add(keras_genomics.layers.RevCompConv1D(filters=100, kernel_size=model.layers[-1].output_shape[1],
        activation="relu", kernel_initializer="he_normal"))

    model.add(keras_genomics.layers.RevCompConv1D(filters=2, kernel_size=model.layers[-1].output_shape[1],
                activation="softmax", kernel_initializer="he_normal"))
    model.add(keras.layers.Flatten())

    return model



def fit_cnn(X, y, codata=None):
    import keras
    import numpy as np

    flank = X.shape[1]//2

    learning_rate = 0.0001
    mini_batch_size = 64
    n_epochs = 1000
    patience = 5
    interval = 1

    if codata is None: 
        codata=np.empty((X.shape[0],0))
    model = create_combo(codata_shape=(codata.shape[1],))

    # 90-10 split into training and validation
    ix = np.random.choice(range(X.shape[0]), round(0.9*X.shape[0]), replace=False)

    X_val, codata_val, y_val = np.delete(X, ix, axis=0), np.delete(codata, ix, axis=0), np.delete(y, ix, axis=0)
    X_train, codata_train, y_train = X[ix], codata[ix], y[ix]

    optim = keras.optimizers.Adam(lr = learning_rate)
    model.compile(optimizer=optim,loss='categorical_crossentropy', metrics=['accuracy'])

    ival = scars_diagnostics.IntervalEvaluation(validation_data=([X_val, codata_val], y_val), flank=flank, patience=patience, interval=interval)
    model.fit(x = [X_train, codata_train], y = y_train, epochs=n_epochs, batch_size = mini_batch_size, callbacks = [ival])

    return model 



def fit_calibration(model, seq, y, weights, flank, codata=None):
    # isotonic regression
    from sklearn.isotonic import IsotonicRegression 
    import numpy as np

    if codata is None:
        codata = np.empty((seq.shape[0],0))

    Z = np.log(model.predict([seq, codata]))

    Z_ref = np.nansum(Z * seq[:, flank, :], axis=1).astype(float)
    Z_ref[np.isinf(Z_ref)] = -100 # isotonic regression does not allow -inf

    y_bin = 1 * (np.argmax(y, axis=1) == np.argmax(seq[:, flank, :], axis=1))

    ir = IsotonicRegression(y_min=0, y_max=1, out_of_bounds='clip')
    ir.fit(Z_ref, y_bin, sample_weight=weights)

    return ir



def prep_data(indices, y, multiplier, seq_in, codata_in=None, expand=False):
    import numpy as np

    if codata_in is None:
        seq = seq_in[indices]
        if expand:
            ix = np.repeat(range(len(indices)), multiplier)
            seq_out = seq[ix]
            y_out = y[ix]
            return seq_out, y_out 
        else:
            ix = np.where(multiplier > 0)
            weight_out = multiplier[ix]
            seq_out = seq[ix]
            y_out = y[ix]
            return seq_out, y_out, weight_out
    else:
        seq = seq_in[indices]
        codata = codata_in[indices]
        if expand:
            ix = np.repeat(range(len(indices)), multiplier)
            seq_out = seq[ix]
            codata_out = codata[ix]
            y_out = y[ix]
            return filter_nan(seq_out, codata_out, y_out)
        else:
            ix = np.where(multiplier > 0)
            weight_out = multiplier[ix]
            seq_out = seq[ix]
            codata_out = codata[ix]
            y_out = y[ix]
            return filter_nan(seq_out, codata_out, y_out, weight_out)



def filter_nan(X, codata, y, weights=None):
    import numpy as np

    ix_without_nan = np.where(~np.any(np.isnan(codata), axis=1))[0]
    X_out, codata_out, y_out = X[ix_without_nan], codata[ix_without_nan], y[ix_without_nan]

    if weights is not None:
        weights_out = weights[ix_without_nan]
        return X_out, codata_out, y_out, weights_out

    return X_out, codata_out, y_out


# function that anchors on central A, C 
# while keras genomics could make this unnecessary by using converted dense layers, 
# it is still necesary to concatenate codata to penultimate layer
def anchor_on_AC (X, y=None):
    import numpy as np

    flank = X.shape[1]//2
    to_rev = (np.sum(X[:, flank, :] * [0,0,1,1], axis=1) == 1)

    X[to_rev] = X[to_rev,::-1,::-1]
    if y is not None:
        y[to_rev] = y[to_rev,::-1]

    return 



