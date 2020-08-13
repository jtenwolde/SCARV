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


def fit_cnn_seq_only(X, y, flank):
    import keras

    learning_rate = 0.0001
    mini_batch_size = 64
    n_epochs = 1000
    patience = 5
    interval = 1

    model = create_genomics_cnn_seqonly()

    # 90-10 split into training and validation
    ix = np.random.choice(range(X.shape[0]), round(0.1*X.shape[0]), replace=False)

    X_train, y_train = np.delete(X, ix, axis=0), np.delete(y, ix, axis=0)
    X_val, y_val = X[ix], y[ix]

    optim = keras.optimizers.Adam(lr = learning_rate)
    model.compile(optimizer=optim,loss='categorical_crossentropy', metrics=['accuracy'])

    ival = scars_diagnostics.IntervalEvaluation(validation_data=(X_val, y_val), flank=flank, patience=patience, interval=interval)
    model.fit(x = X_train, y = y_train, epochs=n_epochs, batch_size = mini_batch_size, callbacks = [ival])

    return model 


def fit_cnn(X, y, flank, codata):
    import keras

    learning_rate = 0.0001
    mini_batch_size = 64
    n_epochs = 1000
    patience = 5
    interval = 1

    model = create_combo(codata_shape=(codata.shape[1],))

    # 90-10 split into training and validation
    ix = np.random.choice(range(X.shape[0]), round(0.1*X.shape[0]), replace=False)

    X_train, codata_train, y_train = np.delete(X, ix, axis=0), np.delete(codata, ix, axis=0), np.delete(y, ix, axis=0)
    X_val, codata_val, y_val = X[ix], codata[ix], y[ix]

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
        Z = np.log(model.predict(seq))
    else:
        Z = np.log(model.predict([seq, codata]))

    Z_ref = np.sum(Z * seq[:, flank, :], axis=1).astype(float)
    y_bin = 1 * (np.argmax(y, axis=1) == np.argmax(seq[:, flank, :], axis=1))

    ir = IsotonicRegression(y_min=0, y_max=1, out_of_bounds='clip')
    ir.fit(Z_ref, y_bin, sample_weight=weights)

    return ir


def filter_nan(X, codata, y):
    import numpy as np

    ix_without_nan = np.where(~np.any(np.isnan(codata), axis=1))[0]
    X_out, codata_out, y_out = X[ix_without_nan], codata[ix_without_nan], y[ix_without_nan]

    return X_out, codata_out, y_out


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
            return [*filter_nan(seq_out, codata_out, y_out), weight_out]



